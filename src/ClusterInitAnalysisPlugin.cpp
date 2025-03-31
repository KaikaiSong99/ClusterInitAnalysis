#include "ClusterInitAnalysisPlugin.h"

#include "PointData/PointData.h"
#include <ClusterData/ClusterData.h>
#include <format>
#include <iostream>
#include <event/Event.h>

#include <QtCore>
#include <QDebug>

Q_PLUGIN_METADATA(IID "nl.BioVault.ClusterInitAnalysisPlugin")

using namespace mv;
using namespace mv::plugin;

// Initialize the random number generator
QRandomGenerator ClusterInitAnalysisPlugin::Point::rng;

// Initialize the dimension names
std::vector<QString> ClusterInitAnalysisPlugin::Point::dimensionNames = std::vector<QString>({"x_pos", "y_pos" , "heading_rad", "x_heading", "y_heading", "vel_position", "vel_heading" });

// Initialize the number of dimensions in the point
std::uint32_t ClusterInitAnalysisPlugin::Point::numberOfDimensions = 7;

// Initialize the maximum velocity
float ClusterInitAnalysisPlugin::Point::maximumVelocity = 1.0f;

ClusterInitAnalysisPlugin::ClusterInitAnalysisPlugin(const PluginFactory* factory) :
    AnalysisPlugin(factory),
    _settingsAction(),
    _points(),
    _pointHeadings()
{
}

void ClusterInitAnalysisPlugin::init()
{
    // Create example output dataset (a points dataset which is derived from the input points dataset) and set the output dataset
    //setOutputDataset(mv::data().createDerivedDataset("Output Data", getInputDataset(), getInputDataset()));
    QString clusterName;
    // Retrieve the input dataset for our specific data type (in our case points)
    // The ManiVault core sets the input dataset reference when the plugin is created
    const auto inputPoints = getInputDataset<Points>();

    int mnist_num_points = 60000;
    // check if the number of points are 44789, if so, it is the cell dataset
    // else check if the number of points are 70000, if so, it is the MNIST dataset
    // otherwise, print an error message
    if (inputPoints->getNumPoints() == 44789) {
        clusterName = "Cell_Cluster";
    }
    else if (inputPoints->getNumPoints() == mnist_num_points) {
        clusterName = "MNIST_Cluster";
    }
    else {
        qDebug() << "Error: Neither Fetal Intestine Immune Cells nor MNIST";
        return;
    }
    setOutputDataset(mv::data().createDataset<Clusters>("Cluster", clusterName, getInputDataset()));
    
    //const auto inputPoints  = getInputDataset<Points>();

    // Retrieve the output dataset for our specific data type (in our case points)
    auto outputPoints = getOutputDataset<Clusters>();

    const auto updatePoints = [this, clusterName, mnist_num_points]() {
        const auto inputPoints = getInputDataset<Points>();
        QVector<QVector<int>> labelsPerClass;
        int numPoints = inputPoints->getNumPoints();

        QStringList clusterNames;
        QVector<QRgb> colors;
        std::string name;

        if (clusterName == "MNIST_Cluster") {
            clusterNames = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9" };
            colors = { qRgb(245, 149, 28), qRgb(166, 180, 166), qRgb(8, 64, 11), qRgb(205, 62, 203), qRgb(198, 116, 108), qRgb(109, 178, 250), qRgb(142, 144, 157), qRgb(29, 126, 24), qRgb(35, 38, 84), qRgb(230, 24, 23) };
            // format 10000 in mnist_10000_labels with variable mnist_num_points
            //name = "D:/TUDelftComputerScience/Thesis/Code/ManiVaultFork/LineViewPlugin/res/data/mnist_" + std::to_string(mnist_num_points) + "_labels.bin";
            name = ":data/mnist_" + std::to_string(mnist_num_points) + "_labels.bin";
            //name = "D:/TUDelftComputerScience/Thesis/Code/ManiVaultFork/LineViewPlugin/res/data/mnist_10000_labels.bin";
        }
            
        else if (clusterName == "Cell_Cluster") {
            clusterNames = { "CD56+CD8a-NK MC17", "CD8a_CD56+NK MC15", "CD56+ILC3 MC6", "ILC2 MC4", "CD161-ILC MC10", "CD56dimCD8a+NK MC16", "atypical cells", "CD56-ILC3MC8 ",
                             "CD8a-MC12", "CD27+MC12", "CD161-ILC3 or LTi MC9", "NKp44+ILC3", "ILC1 MC11", "CD8a+MC12", "unclassified", "CD34+MC1-2", "CD45RA+ILC3 MC7", "CD7_NK MC18" };
            colors = { qRgb(203, 221, 241), qRgb(51, 153, 255), qRgb(8, 64, 11), qRgb(205, 62, 203), qRgb(198, 116, 108), qRgb(109, 178, 250), qRgb(227, 119, 194), qRgb(29, 126, 24), qRgb(230, 24, 23),
                       qRgb(169, 18, 23), qRgb(169, 254, 170), qRgb(116, 240, 160), qRgb(245, 149, 28), qRgb(202, 205, 36), qRgb(7, 48, 106), qRgb(35, 38, 84), qRgb(166, 180, 166), qRgb(142, 144, 157) };
            //name = "D:/TUDelftComputerScience/Thesis/Code/ManiVaultFork/LineViewPlugin/res/data/cell_all_labels.bin";
            name = ":data/cell_all_labels.bin";
        }
            
        

        // initialize the labelPerPoint with the size of numPoints of all 0s
        std::vector<int> labelPerPoint(numPoints, 0);
        // read by creating a Qfile (:/ relative paths works only for QFile, not from std::)
        QFile labelFile(name.c_str());
        if (!labelFile.open(QIODevice::ReadOnly)) {
            qDebug() << "Failed to open input file.";
            return;
        }
        QDataStream in(&labelFile);
        in.readRawData(reinterpret_cast<char*>(labelPerPoint.data()), numPoints * sizeof(int));
        labelFile.close();

        std::vector<int> uniqueLabels = labelPerPoint;
        std::sort(uniqueLabels.begin(), uniqueLabels.end());
        uniqueLabels.erase(std::unique(uniqueLabels.begin(), uniqueLabels.end()), uniqueLabels.end());
        labelsPerClass.resize(uniqueLabels.size());

        for (int i = 0; i < numPoints; ++i) {
            labelsPerClass[labelPerPoint[i]].push_back(labelPerPoint[i]);
        }

        auto outputPoints = getOutputDataset<Clusters>();

        int indicesCount = 0;
        Cluster testCluster;
        for (int j = 0; j < labelsPerClass.size(); j++) {
            Cluster cluster;
            std::vector<std::seed_seq::result_type> indices;
            for (int i = 0; i < labelsPerClass[j].size(); i++)
            {
                indices.push_back(indicesCount + i);
            }
            indicesCount += labelsPerClass[j].size();
            if (labelsPerClass.size() == 1) cluster.setName("All Cells");
            else cluster.setName(clusterNames[j]);
            qDebug() << "cluster name: " << clusterNames[j];
                
            cluster.setColor(colors[j]);
            cluster.setIndices(indices);
            //if (j == 0) testCluster = cluster;
            //else outputPoints->addCluster(cluster);
            outputPoints->addCluster(cluster);
        }
        //outputPoints->addCluster(testCluster);
        
        // Inform the core (and thus others) that the data changed
        events().notifyDatasetDataChanged(outputPoints);
    };

    // Initializes the points
    const auto initializePoints = [this, inputPoints, updatePoints]() {
        // Update the points
        updatePoints();
    };

    // Initialize our points
    initializePoints();

    // Register for points datasets events using a custom callback function
    _eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DatasetAdded));
    _eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DatasetDataChanged));
    _eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DatasetRemoved));
    //_eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DatasetDataSelectionChanged));

    _eventListener.registerDataEventByType(PointType, std::bind(&ClusterInitAnalysisPlugin::onDataEvent, this, std::placeholders::_1));
}

void ClusterInitAnalysisPlugin::onDataEvent(mv::DatasetEvent* dataEvent)
{
    // The data event has a type so that we know what type of data event occurred (e.g. data added, changed, removed, renamed, selection changes)
    switch (dataEvent->getType()) {

        // A points dataset was added
        case EventType::DatasetAdded:
        {
            // Cast the data event to a data added event
            const auto dataAddedEvent = static_cast<DatasetAddedEvent*>(dataEvent);

            // Get the GUI name of the added points dataset and print to the console
            qDebug() << dataAddedEvent->getDataset()->getGuiName() << "was added";

            break;
        }

        // Points dataset data has changed
        case EventType::DatasetDataChanged:
        {
            // Cast the data event to a data changed event
            const auto dataChangedEvent = static_cast<DatasetDataChangedEvent*>(dataEvent);

            // Get the GUI name of the points dataset of which the data changed and print to the console
            qDebug() << dataChangedEvent->getDataset()->getGuiName() << "data changed";

            break;
        }

        // Points dataset data was removed
        case EventType::DatasetRemoved:
        {
            // Cast the data event to a data removed event
            const auto dataRemovedEvent = static_cast<DatasetRemovedEvent*>(dataEvent);

            // Get the GUI name of the removed points dataset and print to the console
            qDebug() << dataRemovedEvent->getDataset()->getGuiName() << "was removed";

            break;
        }

        // Points dataset selection has changed
        case EventType::DatasetDataSelectionChanged:
        {
            // Cast the data event to a data selection changed event
            const auto dataSelectionChangedEvent = static_cast<DatasetDataSelectionChangedEvent*>(dataEvent);

            // Get points dataset
            const auto& changedDataSet = dataSelectionChangedEvent->getDataset();

            // Get the selection set that changed
            const auto selectionSet = changedDataSet->getSelection<Points>();

            // Print to the console
            qDebug() << changedDataSet->getGuiName() << "selection has changed";

            break;
        }

        default:
            break;
    }
}

AnalysisPlugin* ClusterInitAnalysisPluginFactory::produce()
{
    // Return a new instance of the example analysis plugin
    return new ClusterInitAnalysisPlugin(this);
}

mv::DataTypes ClusterInitAnalysisPluginFactory::supportedDataTypes() const
{
    DataTypes supportedTypes;

    // This example analysis plugin is compatible with points datasets
    supportedTypes.append(PointType);

    return supportedTypes;
}

mv::gui::PluginTriggerActions ClusterInitAnalysisPluginFactory::getPluginTriggerActions(const mv::Datasets& datasets) const
{
    PluginTriggerActions pluginTriggerActions;

    const auto getPluginInstance = [this](const Dataset<Points>& dataset) -> ClusterInitAnalysisPlugin* {
        return dynamic_cast<ClusterInitAnalysisPlugin*>(plugins().requestPlugin(getKind(), { dataset }));
    };

    const auto numberOfDatasets = datasets.count();

    if (numberOfDatasets >= 1 && PluginFactory::areAllDatasetsOfTheSameType(datasets, PointType)) {
        auto pluginTriggerAction = new PluginTriggerAction(const_cast<ClusterInitAnalysisPluginFactory*>(this), this, "Cluster Init Analysis", "Perform an Cluster Init Analysis", getIcon(), [this, getPluginInstance, datasets](PluginTriggerAction& pluginTriggerAction) -> void {
            for (auto dataset : datasets)
                getPluginInstance(dataset);
            });

        pluginTriggerActions << pluginTriggerAction;
    }

    return pluginTriggerActions;
}

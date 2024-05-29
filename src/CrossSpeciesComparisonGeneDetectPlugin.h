#pragma once

#include <ViewPlugin.h>

#include <Dataset.h>
#include <widgets/DropWidget.h>

#include <PointData/PointData.h>

#include <QWidget>

/** All plugin related classes are in the ManiVault plugin namespace */
using namespace mv::plugin;

/** Drop widget used in this plugin is located in the ManiVault gui namespace */
using namespace mv::gui;

/** Dataset reference used in this plugin is located in the ManiVault util namespace */
using namespace mv::util;

class QLabel;

class CrossSpeciesComparisonGeneDetectPlugin : public ViewPlugin
{
    Q_OBJECT

public:

    /**
     * Constructor
     * @param factory Pointer to the plugin factory
     */
    CrossSpeciesComparisonGeneDetectPlugin(const PluginFactory* factory);

    /** Destructor */
    ~CrossSpeciesComparisonGeneDetectPlugin() override = default;
    
    /** This function is called by the core after the view plugin has been created */
    void init() override;
    void modifyTableData(QStandardItemModel* model);
    /**
     * Invoked when a data event occurs
     * @param dataEvent Data event which occurred
     */
    void onDataEvent(mv::DatasetEvent* dataEvent);

protected:
    QTableView           *_tableView;                /** Table view for the data */
    //DropWidget*             _dropWidget;                /** Widget for drag and drop behavior */
    //mv::Dataset<Points>   _points;                    /** Points smart pointer */
   // QString                 _currentDatasetName;        /** Name of the current dataset */
    //QLabel*                 _currentDatasetNameLabel;   /** Label that show the current dataset name */
};

/**
 * CrossSpeciesComparisonGeneDetect plugin factory class
 *
 * Note: Factory does not need to be altered (merely responsible for generating new plugins when requested)
 */
class CrossSpeciesComparisonGeneDetectPluginFactory : public ViewPluginFactory
{
    Q_INTERFACES(mv::plugin::ViewPluginFactory mv::plugin::PluginFactory)
    Q_OBJECT
    Q_PLUGIN_METADATA(IID   "studio.manivault.CrossSpeciesComparisonGeneDetectPlugin"
                      FILE  "CrossSpeciesComparisonGeneDetectPlugin.json")

public:

    /** Default constructor */
    CrossSpeciesComparisonGeneDetectPluginFactory() {}

    /** Destructor */
    ~CrossSpeciesComparisonGeneDetectPluginFactory() override {}
    
    /** Creates an instance of the example view plugin */
    ViewPlugin* produce() override;

    /** Returns the data types that are supported by the example view plugin */
    mv::DataTypes supportedDataTypes() const override;

    /**
     * Get plugin trigger actions given \p datasets
     * @param datasets Vector of input datasets
     * @return Vector of plugin trigger actions
     */
    PluginTriggerActions getPluginTriggerActions(const mv::Datasets& datasets) const override;
};

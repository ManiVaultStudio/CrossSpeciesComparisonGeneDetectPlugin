#pragma once
#include <actions/WidgetAction.h>
#include <actions/IntegralAction.h>
#include <actions/OptionAction.h>
#include <actions/OptionsAction.h>
#include <actions/ToggleAction.h>
#include "actions/DatasetPickerAction.h"
#include "PointData/PointData.h"
#include "ClusterData/ClusterData.h"
#include "event/EventListener.h"
#include "actions/Actions.h"
#include "Plugin.h"
#include "DataHierarchyItem.h"
#include "Set.h"
#include <AnalysisPlugin.h>
#include <memory>
#include <algorithm>    
#include <QDebug>
#include <QLabel>
#include <QComboBox>
#include <QGroupBox>
#include <QPushButton>
#include <QGridLayout>
#include <QFormLayout>
#include <QString>
#include <string>
#include <event/Event.h>
#include <QDebug>
#include <QLabel>
#include <string>
#include "actions/VariantAction.h"
#include "actions/GroupAction.h"
using namespace mv::gui;
class QMenu;
class CrossSpeciesComparisonGeneDetectPlugin;

class FetchMetaData;
namespace mv
{
    class CoreInterface;
}

class SettingsAction : public WidgetAction
{
public:
    class OptionSelectionAction : public GroupAction
    {
    protected:
        class Widget : public mv::gui::WidgetActionWidget {
        public:
            Widget(QWidget* parent, OptionSelectionAction* optionSelectionAction);
        };

        QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
            return new OptionSelectionAction::Widget(parent, this);
        };

    public:
        OptionSelectionAction(SettingsAction& SettingsAction);

    protected:
        SettingsAction& _settingsAction;
    };



protected:

    class Widget : public mv::gui::WidgetActionWidget {
    public:
        Widget(QWidget* parent, SettingsAction* SettingsAction);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new SettingsAction::Widget(parent, this);
    };

public:
    SettingsAction(CrossSpeciesComparisonGeneDetectPlugin& CrossSpeciesComparisonGeneDetectPlugins);

public: // Action getters

    VariantAction& getTableModelAction() { return _tableModel; }
    StringAction& getSelectedGeneAction() { return _selectedGene; }
    StringAction&  getSelectedRowIndexAction() { return _selectedRowIndex; }
    DatasetPickerAction& getFilteringTreeDatasetAction() { return _filteringTreeDataset; }
    OptionSelectionAction& getOptionSelectionAction() { return _optionSelectionAction; }
    TriggerAction& getStartComputationTriggerAction() { return _startComputationTriggerAction; }
    DatasetPickerAction& getReferenceTreeDatasetAction() { return _referenceTreeDataset; }
    DatasetPickerAction& getMainPointsDataset() { return _mainPointsDataset; }
    DatasetPickerAction& getEmbeddingDataset() { return _embeddingDataset; }
    DatasetPickerAction& getSpeciesNamesDataset() { return _speciesNamesDataset; }
    DatasetPickerAction& getClusterNamesDataset() { return _clusterNamesDataset; }
    VariantAction& getFilteredGeneNames() { return _filteredGeneNamesVariant; }
    IntegralAction& getTopNGenesFilter() { return _topNGenesFilter; }
    StringAction& getGeneNamesConnection() { return _geneNamesConnection; }
    TriggerAction& getCreateRowMultiSelectTree() { return _createRowMultiSelectTree; }
    ToggleAction& getPerformGeneTableTsneAction() { return _performGeneTableTsneAction; }
    IntegralAction& getTsnePerplexity() { return _tsnePerplexity; }
    OptionsAction& getHiddenShowncolumns() { return _hiddenShowncolumns; }
    OptionAction& getScatterplotColorOption() { return _scatterplotColorOption; }


    void populatePointData(QString& datasetId, std::vector<float>& pointVector, int& numPoints, int& numDimensions, std::vector<QString>& dimensionNames);
    void populateClusterData(QString& datasetId, std::map<QString, std::pair<QColor, std::vector<int>>>& clusterMap);

    double* condensedDistanceMatrix(std::vector<float>& items);
    std::string mergeToNewick(int* merge, int numOfLeaves);
    QString createJsonTreeFromNewick(QString tree, std::vector<QString> leafNames);
private:
    QVariant createModelFromData(const QStringList& returnGeneList, const std::map<QString, std::map<QString, float>>& map, const QString& treeDatasetId, const float& treeSimilarityScore, const std::map<QString, std::vector<QString>>& geneCounter, const int& n);
    QVariant findTopNGenesPerCluster(const std::map<QString, std::map<QString, float>>& map, int n, QString datasetId, float treeSimilarityScore);

public: // Serialization

    /**
     * Load widget action from variant map
     * @param Variant map representation of the widget action
     */
    void fromVariantMap(const QVariantMap& variantMap) override;

    /**
     * Save widget action to variant map
     * @return Variant map representation of the widget action
     */
    QVariantMap toVariantMap() const override;

protected:
    CrossSpeciesComparisonGeneDetectPlugin& _crossSpeciesComparisonGeneDetectPlugin;
    VariantAction                 _tableModel;
    StringAction                  _selectedGene;
    DatasetPickerAction          _filteringTreeDataset;
    StringAction                _selectedRowIndex;
    OptionSelectionAction         _optionSelectionAction;
    TriggerAction              _startComputationTriggerAction;
    DatasetPickerAction    _referenceTreeDataset;
    std::map<QString, std::map<QString, float>> _clusterGeneMeanExpressionMap;
    DatasetPickerAction    _mainPointsDataset;
    DatasetPickerAction    _speciesNamesDataset;
    DatasetPickerAction    _clusterNamesDataset;
    DatasetPickerAction    _embeddingDataset;
    std::map<QString, std::map<QString, float>> _clusterNameToGeneNameToExpressionValue;
    VariantAction           _filteredGeneNamesVariant;
    IntegralAction          _topNGenesFilter;
    StringAction           _geneNamesConnection;
    TriggerAction         _createRowMultiSelectTree;
    ToggleAction            _performGeneTableTsneAction;
    IntegralAction         _tsnePerplexity;
    OptionsAction          _hiddenShowncolumns;
    Dataset<Points>        _tsneDataset;
    Dataset<Points>        _selectedPointsDataset;
    Dataset<Clusters>        _tsneDatasetSpeciesColors;
    Dataset<Clusters>        _tsneDatasetClusterColors;
    OptionAction           _scatterplotColorOption;
};
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
#include "QStatusBar"
#include <widgets/FlowLayout.h>
#include <QSplitter>

using namespace mv::gui;
class QMenu;
class CrossSpeciesComparisonGeneDetectPlugin;

class FetchMetaData;
namespace mv
{
    class CoreInterface;
}

struct Statistics {
    float meanSelected;
    float medianSelected;
    float modeSelected;
    float rangeSelected;
    int countSelected;
    float meanNonSelected;
    float medianNonSelected;
    float modeNonSelected;
    float rangeNonSelected;
    int countNonSelected;

};
struct SpeciesColorCountStorage {
    QColor color;
    int selectedCellsCount;
    int nonSelectedCellsCount;
};

struct StatisticsSingle {
    float meanVal;
    float medianVal;
    float modeVal;
    float rangeVal;
    int countVal;

};

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
    DatasetPickerAction& getFilteringEditTreeDatasetAction() { return _filteringEditTreeDataset; }
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
    DatasetPickerAction& getScatterplotEmbeddingPointsUMAPOption() { return _scatterplotEmbeddingPointsUMAPOption; }
    OptionAction& getScatterplotReembedColorOption() { return _scatterplotReembedColorOption; }
    StringAction& getSelctedSpeciesVals() { return _selectedSpeciesVals; }
    TriggerAction& getRemoveRowSelection() { return _removeRowSelection; }
    StringAction& getStatusColorAction() { return _statusColorAction; }
    OptionAction& getTypeofTopNGenes() { return _typeofTopNGenes; }
    ToggleAction& getUsePreComputedTSNE() { return _usePreComputedTSNE; }
    //tsne relatedDatasets
    /*
        Dataset<Points>        _selectedPointsTSNEDataset;
    Dataset<Points>        _selectedPointsDataset;
    Dataset<Points>        _selectedPointsEmbeddingDataset;

    Dataset<Clusters>        _tsneDatasetSpeciesColors;
    Dataset<Clusters>        _tsneDatasetClusterColors;
    Dataset<Points>        _tsneDatasetExpressionColors;
    */

    Dataset<Points>& getSelectedPointsTSNEDataset() { return _selectedPointsTSNEDataset; }
    Dataset<Points>& getSelectedPointsDataset() { return _selectedPointsDataset; }
    Dataset<Points>& getSelectedPointsEmbeddingDataset() { return _selectedPointsEmbeddingDataset; }

    Dataset<Clusters>& getTsneDatasetSpeciesColors() { return _tsneDatasetSpeciesColors; }
    Dataset<Clusters>& getTsneDatasetClusterColors() { return _tsneDatasetClusterColors; }
    Dataset<Points>& getTsneDatasetExpressionColors() { return _tsneDatasetExpressionColors; }
    std::vector<std::seed_seq::result_type>& getSelectedIndicesFromStorage() { return _selectedIndicesFromStorage; }
    Dataset<Points> & getFilteredUMAPDatasetPoints() { return _filteredUMAPDatasetPoints; }
    Dataset<Points> & getFilteredUMAPDatasetColors() { return _filteredUMAPDatasetColors; }
    QStatusBar* getStatusBarActionWidget() const { return _statusBarActionWidget; }
    QStringList& getInitColumnNames() { return _initColumnNames; }
    mv::gui::FlowLayout* getSelectedCellClusterInfoStatusBar() const { return _selectedCellClusterInfoStatusBar; }
    QTableView* getTableView() const { return _tableView; }
    QTableView* getSelectionDetailsTable() const { return _selectionDetailsTable; }
    std::map<QString, SpeciesColorCountStorage> & getSelectedSpeciesCellCountMap() { return _selectedSpeciesCellCountMap; }
    QHBoxLayout* getTableSplitter() const { return _splitter; }

    void computeGeneMeanExpressionMap();
    void populatePointDataConcurrently(QString datasetId, const std::vector<float>& pointVector, int numPoints, int numDimensions, std::vector<QString> dimensionNames);
    void populatePointData(QString& datasetId, std::vector<float>& pointVector, int& numPoints, int& numDimensions, std::vector<QString>& dimensionNames);
    void populateClusterData(QString& datasetId, std::map<QString, std::pair<QColor, std::vector<int>>>& clusterMap);

    double* condensedDistanceMatrix(const std::vector<float>& items);
    std::string mergeToNewick(int* merge, int numOfLeaves);
    QString createJsonTreeFromNewick(QString tree, std::vector<QString> leafNames, std::map <QString, Statistics> speciesMeanValues);
private:
    QVariant createModelFromData(const QSet<QString>& returnGeneList, const std::map<QString, std::map<QString, Statistics>>& map, const QString& treeDatasetId, const float& treeSimilarityScore, const std::map<QString, std::vector<QString>>& geneCounter, const int& n);
    QVariant findTopNGenesPerCluster(const std::map<QString, std::map<QString, Statistics>>& map, int n, QString datasetId, float treeSimilarityScore);
    void updateSelectedSpeciesCounts(QJsonObject& node, const std::map<QString, int>& speciesCountMap);
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
    DatasetPickerAction          _filteringEditTreeDataset;
    StringAction                _selectedRowIndex;
    OptionSelectionAction         _optionSelectionAction;
    TriggerAction              _startComputationTriggerAction;
    DatasetPickerAction    _referenceTreeDataset;
    std::map<QString, std::map<QString, float>> _clusterGeneMeanExpressionMap;
    DatasetPickerAction    _mainPointsDataset;
    DatasetPickerAction    _speciesNamesDataset;
    DatasetPickerAction    _clusterNamesDataset;
    DatasetPickerAction    _embeddingDataset;
    std::map<QString, std::map<QString, Statistics>> _clusterNameToGeneNameToExpressionValue;
    VariantAction           _filteredGeneNamesVariant;
    IntegralAction          _topNGenesFilter;
    StringAction           _geneNamesConnection;
    TriggerAction         _createRowMultiSelectTree;
    ToggleAction            _performGeneTableTsneAction;
    IntegralAction         _tsnePerplexity;
    OptionsAction          _hiddenShowncolumns;
    StringAction            _selectedSpeciesVals;
    OptionAction                    _typeofTopNGenes;

    Dataset<Points>        _selectedPointsTSNEDataset;
    Dataset<Points>        _selectedPointsDataset;
    Dataset<Points>        _selectedPointsEmbeddingDataset;
    Dataset<Points>        _filteredUMAPDatasetPoints;
    Dataset<Points>        _filteredUMAPDatasetColors;

    Dataset<Clusters>        _tsneDatasetSpeciesColors;
    Dataset<Clusters>        _tsneDatasetClusterColors;
    Dataset<Points>        _tsneDatasetExpressionColors;
    TriggerAction          _removeRowSelection;
    DatasetPickerAction           _scatterplotEmbeddingPointsUMAPOption;
    OptionAction           _scatterplotReembedColorOption;
    StringAction    _statusColorAction;
    std::vector<std::seed_seq::result_type> _selectedIndicesFromStorage;
    QStatusBar*                     _statusBarActionWidget;
    mv::gui::FlowLayout*            _selectedCellClusterInfoStatusBar;
    //mv::gui::FlowLayout     _clustersLayout;
    QStringList _initColumnNames;
    ToggleAction                  _usePreComputedTSNE;
    QLabel* _currentCellSelectionClusterInfoLabel;
    std::map<QString, SpeciesColorCountStorage>       _selectedSpeciesCellCountMap;

    QTableView* _tableView;                /** Table view for the data */
    QTableView* _selectionDetailsTable;    /** Table view for the selection details */
    QHBoxLayout* _splitter;

};
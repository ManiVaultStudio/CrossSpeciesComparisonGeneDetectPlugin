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
#include <QTableView>
#include <QListView>
#include <vector>
#include <set>
#include <QLineEdit>
#include <QColor>
using namespace mv::gui;
class QMenu;
class CrossSpeciesComparisonGeneDetectPlugin;

class FetchMetaData;
namespace mv
{
    class CoreInterface;
}

class CustomLineEdit : public QLineEdit {
    Q_OBJECT

public:
    explicit CustomLineEdit(QWidget* parent = nullptr) : QLineEdit(parent) {
        QIcon defocusIcon = Application::getIconFont("FontAwesome").getIcon("times-circle");
        _action = addAction(defocusIcon, QLineEdit::TrailingPosition);
        _action->setVisible(false); // Initially hide the action icon

        connect(_action, &QAction::triggered, this, &CustomLineEdit::defocusLineEdit);
    }

signals:
    void textboxSelectedForTyping();
    void textboxDeselectedNotTypingAnymore();
    //void textValueChanged();

protected:
    void focusInEvent(QFocusEvent* e) override {
        _action->setVisible(true);
        emit textboxSelectedForTyping();
        QLineEdit::focusInEvent(e);
    }

    void focusOutEvent(QFocusEvent* e) override {
        // Only hide the defocus icon and emit the deselection signal if the text box is empty
        if (this->text().isEmpty()) {
            _action->setVisible(false);
            emit textboxDeselectedNotTypingAnymore();
        }
        QLineEdit::focusOutEvent(e);
    }

 


private slots:
    void defocusLineEdit() {
        this->blockSignals(true); // Block signals to prevent textChanged from being emitted
        this->clear();
        this->blockSignals(false); // Unblock signals

        this->clearFocus(); // Always clear focus when the defocus button is clicked
        // Optionally, clear the text here if needed
        //this->clear();
        _action->setVisible(false); // Hide the defocus icon
        emit textboxDeselectedNotTypingAnymore(); // Emit the deselection signal
    }

private:
    QAction* _action; // Action for the defocus icon
};




struct ClusterOrderContainer {
    int count;
    QColor color;
    QString name;
};


struct SpeciesDetailsStats {
    int rank;
    float meanSelected;
    int countSelected;
    float meanNonSelected;
    int countNonSelected;
    //float meanAll;
    //int countAll;

};

struct Stats {
    
    float meanSelected;
    int countSelected;
    float meanNonSelected;
    int countNonSelected;
    QColor color;
    //float meanAll;
    //int countAll;

};
struct SpeciesColorCountStorageVals {
    QColor color;
    int selectedCellsCount;
    int nonSelectedCellsCount;
    //int allCellsCount;
};

struct StatisticsSingle {
    float meanVal;
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

    VariantAction& getListModelAction() { return _listModel; }
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
    OptionsAction& getSpeciesExplorerInMap() { return _speciesExplorerInMap; }
    TriggerAction& getSpeciesExplorerInMapTrigger() { return _speciesExplorerInMapTrigger; }
    TriggerAction& getRevertRowSelectionChangesToInitial() { return _revertRowSelectionChangesToInitial; }
    ToggleAction& getApplyLogTransformation() { return _applyLogTransformation; }
    OptionAction& getClusterCountSortingType() { return _clusterCountSortingType; }
    IntegralAction& getPerformGeneTableTsnePerplexity() { return _performGeneTableTsnePerplexity; }
    //IntegralAction& setPerformGeneTableTsnePerplexity() { return _performGeneTableTsnePerplexity; }
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
    QTableView* getGeneTableView() const { return _geneTableView; }
    CustomLineEdit* getSearchBox() const { return _searchBox; }
    QTableView* getSelectionDetailsTable() const { return _selectionDetailsTable; }
    std::map<QString, SpeciesColorCountStorageVals> & getSelectedSpeciesCellCountMap() { return _selectedSpeciesCellCountMap; }
    QHBoxLayout* getTableSplitter() const { return _splitter; }
    std::vector<QString>& getCustomOrderClustersFromHierarchy() { return _customOrderClustersFromHierarchy; }
    std::map<QString, std::map<QString, Stats>>& getClusterNameToGeneNameToExpressionValue() { return _clusterNameToGeneNameToExpressionValue; }
    QSet<QString>& getUniqueReturnGeneList() { return _uniqueReturnGeneList; }
    std::vector<QString>& getTotalGeneList() { return _totalGeneList; }
    Dataset<Points>& getGeneSimilarityPoints() { return _geneSimilarityPoints; }
    //std::vector<QString>& getGeneSimilarityClusters() { return _geneSimilarityClusters; }
    Dataset<Clusters>& getGeneSimilarityClusterColoring() { return _geneSimilarityClusterColoring; }

    //bool getErroredOutFlag() const { return _erroredOutFlag; }
    //bool setErrorOutFlag(bool flag) { return _erroredOutFlag = flag; }

    void computeGeneMeanExpressionMap();
    void populatePointDataConcurrently(QString datasetId, const std::vector<float>& pointVector, int numPoints, int numDimensions, std::vector<QString> dimensionNames);
    void populatePointData(QString& datasetId, std::vector<float>& pointVector, int& numPoints, int& numDimensions, std::vector<QString>& dimensionNames);
    void populateClusterData(QString& datasetId, std::map<QString, std::pair<QColor, std::vector<int>>>& clusterMap);
    QStringList getSystemModeColor();
    double* condensedDistanceMatrix(const std::vector<float>& items);
    std::string mergeToNewick(int* merge, int numOfLeaves);
    QString createJsonTreeFromNewick(QString tree, std::vector<QString> leafNames, std::map <QString, Stats> speciesMeanValues);
    void disableActions();
    void enableActions();
    void removeSelectionTableRows(QStringList* selectedLeaves);
    void enableDisableButtonsAutomatically();

    QVariant createModelFromData(const std::map<QString, std::map<QString, Stats>>& map, const std::map<QString, std::vector<QString>>& geneCounter, const std::map<QString, std::vector<std::pair<QString, int>>>& rankingMap,const int& n);
    void findTopNGenesPerCluster();
private:

    void updateSelectedSpeciesCounts(QJsonObject& node, const std::map<QString, int>& speciesCountMap);
    void updateButtonTriggered();
    void setModifiedTriggeredData(QVariant geneListTable);
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
    VariantAction                 _listModel;
    StringAction                  _selectedGene;
    DatasetPickerAction          _filteringEditTreeDataset;
    StringAction                _selectedRowIndex;
    OptionSelectionAction         _optionSelectionAction;
    TriggerAction              _startComputationTriggerAction;
    DatasetPickerAction    _referenceTreeDataset;
    std::map<QString, std::map<QString, std::pair<int,float>>> _clusterGeneMeanExpressionMap;
    DatasetPickerAction    _mainPointsDataset;
    DatasetPickerAction    _speciesNamesDataset;
    DatasetPickerAction    _clusterNamesDataset;
    DatasetPickerAction    _embeddingDataset;
    std::map<QString, std::map<QString, Stats>> _clusterNameToGeneNameToExpressionValue;
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

    Dataset<Points>             _geneSimilarityPoints;
    //std::vector<QString>        _geneSimilarityClusters;
    Dataset<Clusters>           _geneSimilarityClusterColoring;

    TriggerAction          _removeRowSelection;
    TriggerAction           _revertRowSelectionChangesToInitial;
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
    std::map<QString, SpeciesColorCountStorageVals>       _selectedSpeciesCellCountMap;
    OptionsAction                               _speciesExplorerInMap;
    TriggerAction                               _speciesExplorerInMapTrigger;
    QTableView* _geneTableView;                /** Table view for the data */
    QTableView* _selectionDetailsTable;    /** Table view for the selection details */
    QHBoxLayout* _splitter;
    CustomLineEdit* _searchBox;
    ToggleAction    _applyLogTransformation;
    //bool _erroredOutFlag;
    bool _meanMapComputed;
    OptionAction                _clusterCountSortingType;

    std::vector<QString> _customOrderClustersFromHierarchy;
    std::unordered_map<QString, int> _customOrderClustersFromHierarchyMap;
    std::vector<QString>     _totalGeneList;
    QSet<QString>               _uniqueReturnGeneList;
    IntegralAction                _performGeneTableTsnePerplexity;
};
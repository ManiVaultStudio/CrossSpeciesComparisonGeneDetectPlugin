#include "SettingsAction.h"
#include "CrossSpeciesComparisonGeneDetectPlugin.h"
#include<string>  
#include <QFileDialog>
#include <QPageLayout>
#include <QWebEngineView>
#include <CrossSpeciesComparisonTreeData.h>
#include <numeric>   // for std::reduce
#include <execution> // for std::execution::par
//#include "lib/Distance/annoylib.h"
//#include "lib/Distance/kissrandom.h"
#include <QtConcurrent>
//#include "lib/JSONnlohmann/json.hpp"
//#include "lib/Clustering/fastcluster.h"
//#include "lib/NewickComparator/newick_comparator.h" //https://github.com/MaciejSurowiec/Maximum_agreement_subtree_problem
#include <sstream>
#include <stack>
#include <algorithm> // for std::find
#include <mutex>
#include <unordered_set>
#include <iostream>
#include <map>
#include <string>
#include <chrono>
#include <cmath> // Include for std::log
#include <cmath> // Include for std::isnan and std::isinf
#include <limits>
#include <vector>

#include <unordered_map>
#include <QApplication>
#include <QPalette>

using namespace mv;
using namespace mv::gui;





std::map<std::string, std::chrono::high_resolution_clock::time_point> timers;



// Function to sort by count (descending)
bool sortByCount(const ClusterOrderContainer& a, const ClusterOrderContainer& b) {
    return a.count > b.count;
}

// Function to sort by name (ascending)
bool sortByName(const ClusterOrderContainer& a, const ClusterOrderContainer& b) {
    return a.name < b.name;
}

// Prepare a map for custom list sorting
std::unordered_map<QString, int> prepareCustomSortMap(const std::vector<QString>& customOrder) {
    std::unordered_map<QString, int> sortOrderMap;
    for (size_t i = 0; i < customOrder.size(); ++i) {
        sortOrderMap[customOrder[i]] = i;
    }
    return sortOrderMap;
}

// Function to sort by custom list
bool sortByCustomList(const ClusterOrderContainer& a, const ClusterOrderContainer& b,
    const std::unordered_map<QString, int>& sortOrderMap) {
    auto posA = sortOrderMap.find(a.name);
    auto posB = sortOrderMap.find(b.name);
    int indexA = (posA != sortOrderMap.end()) ? posA->second : INT_MAX;
    int indexB = (posB != sortOrderMap.end()) ? posB->second : INT_MAX;
    return indexA < indexB;
}




Stats combineStatisticsSingle(const StatisticsSingle& selected, const StatisticsSingle& nonSelected/*, const StatisticsSingle& allSelected*/) {
    Stats combinedStats;
    combinedStats.meanSelected = selected.meanVal;
    combinedStats.countSelected = selected.countVal;

    //combinedStats.meanAll = allSelected.meanVal;
    //combinedStats.countAll = allSelected.countVal;

    combinedStats.meanNonSelected = nonSelected.meanVal; 
    combinedStats.countNonSelected = nonSelected.countVal; 

    return combinedStats;
}

StatisticsSingle calculateStatistics(const std::vector<float>& data) {
    const int count = data.size();
    if (count == 0) {
        return { 0.0f, 0 }; // Return early if data is empty to avoid division by zero
    }

    float sum = 0.0;
#ifdef _WIN32
    sum = std::reduce(std::execution::par, data.begin(), data.end(), 0.0f);
#else
    sum = std::reduce(data.begin(), data.end(), 0.0f);
#endif

    float mean = std::round((sum / count) * 100.0f) / 100.0f;

    return { mean, count };
}
 float calculateMean(const std::vector<float>&v) {
    if (v.empty())
         return 0.0f;
    
        
        #ifdef _WIN32
         float sum = std::reduce(std::execution::par, v.begin(), v.end(), 0.0f);
    #else
         float sum = std::reduce(v.begin(), v.end(), 0.0f);
    #endif
        
        
        return sum / static_cast<float>(v.size());
    
}


void startCodeTimer(const std::string& message) {
    timers[message] = std::chrono::high_resolution_clock::now();
    std::cout << "Timer started for: " << message << std::endl;
}

void stopCodeTimer(const std::string& message) {
    auto search = timers.find(message);
    if (search != timers.end()) {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - search->second;
        std::cout << "Timer stopped for: " << message << ". Elapsed time: " << elapsed.count() << "s\n";
        timers.erase(search); // Remove the timer from the map
    }
    else {
        std::cout << "Timer not found for: " << message << ". Ensure the timer was started with the exact same message string." << std::endl;
    }
}


/*
float calculateMeanLogTransformed(const std::vector<float>& v) {
    if (v.empty())
        return 0.0f;

    // Use a lambda to filter and log-transform positive values only
    auto logTransform = [](float value) -> float {
        return value > 0.0f ? std::log(value) : 0.0f;
        };

    // Count how many positive values are present
    size_t positiveCount = std::count_if(v.begin(), v.end(), [](float value) { return value > 0.0f; });

    // If there are no positive values, return 0 to avoid division by zero
    if (positiveCount == 0)
        return 0.0f;

    

    #ifdef _WIN32
    float sum = std::transform_reduce(std::execution::par, v.begin(), v.end(), 0.0f, std::plus<>(), logTransform);
    #else
    float sum = std::transform_reduce(v.begin(), v.end(), 0.0f, std::plus<>(), logTransform);
    #endif

    return sum / static_cast<float>(positiveCount);
}
*/
/*
std::string jsonToNewick(const nlohmann::json& node, const std::vector<QString>& species) {
    std::string newick;
    if (node.contains("children")) {
        newick += "(";
        for (const auto& child : node["children"]) {
            newick += jsonToNewick(child, species) + ",";
        }
        newick = newick.substr(0, newick.size() - 1);  // Remove trailing comma
        newick += ")";
    }
    if (node.contains("name")) {
        std::string nodeName = node["name"].get<std::string>();
        auto it = std::find_if(species.begin(), species.end(), [&nodeName](const QString& str) {
            return str.compare(QString::fromStdString(nodeName)) == 0;
            });
        if (it != species.end()) {
            newick += std::to_string(std::distance(species.begin(), it) + 1);  // Indices start from 1
        }
    }
    return newick;
}
*/
bool areSameIgnoreOrder(const QStringList& list1, const QStringList& list2) {
    if (list1.size() != list2.size()) {
        return false;
    }

    QHash<QString, int> list1Counts;
    for (const auto& item : list1) {
        ++list1Counts[item];
    }

    for (const auto& item : list2) {
        if (!list1Counts.contains(item) || --list1Counts[item] < 0) {
            return false;
        }
    }

    return true;
}



int findIndex(const std::vector<std::seed_seq::result_type>& vec, int value) {
    auto it = std::find(vec.begin(), vec.end(), value);
    return (it != vec.end()) ? std::distance(vec.begin(), it) : -1;
}


SettingsAction::SettingsAction(CrossSpeciesComparisonGeneDetectPlugin& CrossSpeciesComparisonGeneDetectPlugin) :
    WidgetAction(&CrossSpeciesComparisonGeneDetectPlugin, "CrossSpeciesComparisonGeneDetectPlugin Settings"),
    _crossSpeciesComparisonGeneDetectPlugin(CrossSpeciesComparisonGeneDetectPlugin),
    _listModel(this, "List Model"),
    _selectedGene(this, "Selected Gene"),
    _filteringEditTreeDataset(this, "Filtering Tree Dataset"),
    _selectedRowIndex(this, "Selected Row Index"),
    _optionSelectionAction(*this),
    _startComputationTriggerAction(this, "Update"),
    _referenceTreeDataset(this, "Reference Tree Dataset"),
    _mainPointsDataset(this, "Main Points Dataset"),
    _embeddingDataset(this, "Embedding Dataset"),
    //_hierarchyBottomClusterDataset(this, "Hierarchy Bottom Cluster Dataset"),
    //_hierarchyMiddleClusterDataset(this, "Hierarchy Middle Cluster Dataset"),
    //_hierarchyTopClusterDataset(this, "Hierarchy Top Cluster Dataset"),
    _speciesNamesDataset(this, "Species Names"),
    _clusterNamesDataset(this, "Cluster Names"),
    //_calculationReferenceCluster(this, "Calculation Reference Cluster"),
    _filteredGeneNamesVariant(this, "Filtered Gene Names"),
    _topNGenesFilter(this, "Top N"),
    _geneNamesConnection(this, "Gene Names Connection"),
    _createRowMultiSelectTree(this, "Create Row MultiSelect Tree"),
    _performGeneTableTsneAction(this, "Perform Gene Table TSNE"),
    _tsnePerplexity(this, "TSNE Perplexity"),
    _hiddenShowncolumns(this, "Hidden Shown Columns"),
    _scatterplotReembedColorOption(this, "Embed Color"),
    _scatterplotEmbeddingPointsUMAPOption(this, "Embedding UMAP Points"),
    _selectedSpeciesVals(this, "Selected Species Vals"),
    _removeRowSelection(this, "DeSelect"),
    _revertRowSelectionChangesToInitial(this, "Revert"),
    _statusColorAction(this, "Status color"),
    _typeofTopNGenes(this, "N Type"),
    _usePreComputedTSNE(this, "Use Precomputed TSNE"),
    _speciesExplorerInMap(this, "Leaves Explorer Options"),
    _speciesExplorerInMapTrigger(this, "Explore"),
    _applyLogTransformation(this, "Gene mapping log"),
    _clusterCountSortingType(this, "Cluster Count Sorting Type"),
    _currentCellSelectionClusterInfoLabel(nullptr)
{
    
    setSerializationName("CSCGDV:CrossSpeciesComparison Gene Detect Plugin Settings");
    _statusBarActionWidget  = new QStatusBar();
    _searchBox = new QLineEdit();
    QIcon searchIcon = Application::getIconFont("FontAwesome").getIcon("search");
    QAction* searchAction = new QAction(_searchBox);
    searchAction->setIcon(searchIcon);
    _searchBox->addAction(searchAction, QLineEdit::LeadingPosition);
    _searchBox->setPlaceholderText("Search ID...");
    _searchBox->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    _searchBox->setMaximumHeight(22);
    //_searchBox->setMinimumWidth(100);
    _searchBox->setMaximumWidth(800);
    _searchBox->setAutoFillBackground(true);
    _searchBox->setStyleSheet("QLineEdit { background-color: white; }");
    _searchBox->setClearButtonEnabled(true);
    _searchBox->setFocusPolicy(Qt::StrongFocus);
    _meanMapComputed = false;
    _statusBarActionWidget->setStatusTip("Status");
    _statusBarActionWidget->setMaximumHeight(22);
    //_statusBarActionWidget->setFixedWidth(120);
    //_statusBarActionWidget->setMinimumWidth(100);
    _statusBarActionWidget->setMaximumWidth(800);
    _statusBarActionWidget->setAutoFillBackground(true);
    _statusBarActionWidget->setSizeGripEnabled(false);

    _geneTableView = new QTableView();
    _selectionDetailsTable = new QTableView();
    _splitter = new QHBoxLayout();
    _geneTableView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    _geneTableView->setSelectionBehavior(QAbstractItemView::SelectRows);
    _geneTableView->setSelectionMode(QAbstractItemView::SingleSelection);
    _geneTableView->setEditTriggers(QAbstractItemView::NoEditTriggers);
    _geneTableView->setAlternatingRowColors(true);
    _geneTableView->setSortingEnabled(true);
    _geneTableView->setShowGrid(true);
    _geneTableView->setGridStyle(Qt::SolidLine);
    _geneTableView->setHorizontalScrollMode(QAbstractItemView::ScrollPerPixel);
    _geneTableView->setVerticalScrollMode(QAbstractItemView::ScrollPerPixel);
    _geneTableView->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _geneTableView->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _geneTableView->setCornerButtonEnabled(false);
    _geneTableView->setWordWrap(false);
    _geneTableView->setTabKeyNavigation(false);
    _geneTableView->setAcceptDrops(false);
    _geneTableView->setDropIndicatorShown(false);
    _geneTableView->setDragEnabled(false);
    _geneTableView->setDragDropMode(QAbstractItemView::NoDragDrop);
    _geneTableView->setDragDropOverwriteMode(false);
    _geneTableView->setAutoScroll(false);
    _geneTableView->setAutoScrollMargin(16);
    _geneTableView->setAutoFillBackground(true);
    _geneTableView->setFrameShape(QFrame::NoFrame);
    _geneTableView->setFrameShadow(QFrame::Plain);
    _geneTableView->setLineWidth(0);
    _geneTableView->setMidLineWidth(0);
    _geneTableView->setFocusPolicy(Qt::NoFocus);
    _geneTableView->setContextMenuPolicy(Qt::NoContextMenu);
    _geneTableView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    _geneTableView->setMinimumSize(QSize(0, 0));
    _geneTableView->setMaximumSize(QSize(16777215, 16777215));
    _geneTableView->setBaseSize(QSize(0, 0));
    _geneTableView->setFocusPolicy(Qt::StrongFocus);
    _geneTableView->setVerticalScrollMode(QAbstractItemView::ScrollPerPixel);

    //only highlight multiple rows if shiuft is pressed
    _geneTableView->setSelectionBehavior(QAbstractItemView::SelectRows);



    _selectionDetailsTable->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    _selectionDetailsTable->setSelectionBehavior(QAbstractItemView::SelectRows);
    _selectionDetailsTable->setSelectionMode(QAbstractItemView::SingleSelection);
    _selectionDetailsTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    _selectionDetailsTable->setAlternatingRowColors(false);
    _selectionDetailsTable->setSortingEnabled(true);
    _selectionDetailsTable->setShowGrid(true);
    _selectionDetailsTable->setGridStyle(Qt::SolidLine);
    _selectionDetailsTable->setHorizontalScrollMode(QAbstractItemView::ScrollPerPixel);
    _selectionDetailsTable->setVerticalScrollMode(QAbstractItemView::ScrollPerPixel);
    _selectionDetailsTable->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _selectionDetailsTable->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _selectionDetailsTable->setCornerButtonEnabled(false);
    _selectionDetailsTable->setWordWrap(false);
    _selectionDetailsTable->setTabKeyNavigation(false);
    _selectionDetailsTable->setAcceptDrops(false);
    _selectionDetailsTable->setDropIndicatorShown(false);
    _selectionDetailsTable->setDragEnabled(false);
    _selectionDetailsTable->setDragDropMode(QAbstractItemView::NoDragDrop);
    _selectionDetailsTable->setDragDropOverwriteMode(false);
    _selectionDetailsTable->setAutoScroll(false);
    _selectionDetailsTable->setAutoScrollMargin(16);
    _selectionDetailsTable->setAutoFillBackground(true);
    _selectionDetailsTable->setFrameShape(QFrame::NoFrame);
    _selectionDetailsTable->setFrameShadow(QFrame::Plain);
    _selectionDetailsTable->setLineWidth(0);
    _selectionDetailsTable->setMidLineWidth(0);
    _selectionDetailsTable->setFocusPolicy(Qt::NoFocus);
    _selectionDetailsTable->setContextMenuPolicy(Qt::NoContextMenu);
    _selectionDetailsTable->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    _selectionDetailsTable->setMinimumSize(QSize(0, 0));
    _selectionDetailsTable->setMaximumSize(QSize(16777215, 16777215));
    _selectionDetailsTable->setBaseSize(QSize(0, 0));
    _selectionDetailsTable->setFocusPolicy(Qt::StrongFocus);
    _selectionDetailsTable->setVerticalScrollMode(QAbstractItemView::ScrollPerPixel);
    _selectionDetailsTable->setSelectionBehavior(QAbstractItemView::SelectRows);







    _selectedCellClusterInfoStatusBar = new mv::gui::FlowLayout();


    _listModel.setSerializationName("CSCGDV:List Model");
    _selectedGene.setSerializationName("CSCGDV:Selected Gene");
    _mainPointsDataset  .setSerializationName("CSCGDV:Main Points Dataset");
    _embeddingDataset.setSerializationName("CSCGDV:Embedding Dataset");
    _speciesNamesDataset.setSerializationName("CSCGDV:Species Names Dataset");
    _clusterNamesDataset.setSerializationName("CSCGDV:Cluster Names Dataset");
    _filteredGeneNamesVariant.setSerializationName("CSCGDV:Filtered Gene Names");
    _topNGenesFilter.setSerializationName("CSCGDV:Top N Genes Filter");
    _filteringEditTreeDataset.setSerializationName("CSCGDV:Filtering Tree Dataset");
    _referenceTreeDataset.setSerializationName("CSCGDV:Reference Tree Dataset");
    _selectedRowIndex.setSerializationName("CSCGDV:Selected Row Index");
    _geneNamesConnection.setSerializationName("CSCGDV:Gene Names Connection");
    _selectedSpeciesVals.setSerializationName("CSCGDV:Selected Species Vals");
    _removeRowSelection.setSerializationName("CSCGDV:Remove Row Selection");
    _removeRowSelection.setDisabled(true);
    _revertRowSelectionChangesToInitial.setSerializationName("CSCGDV:Revert Row Selection Changes To Initial");
    _revertRowSelectionChangesToInitial.setDisabled(true);
    _speciesExplorerInMapTrigger.setSerializationName("CSCGDV:Species Explorer In Map Trigger");
    _speciesExplorerInMapTrigger.setDisabled(true);
    _statusColorAction.setSerializationName("CSCGDV:Status Color");
    _selectedGene.setDisabled(true);
    _selectedGene.setString("");
    _createRowMultiSelectTree.setSerializationName("CSCGDV:Create Row MultiSelect Tree");
    _performGeneTableTsneAction.setSerializationName("CSCGDV:Perform Gene Table TSNE");
    _tsnePerplexity.setSerializationName("CSCGDV:TSNE Perplexity");
    _tsnePerplexity.setMinimum(1);
    _tsnePerplexity.setMaximum(50);
    _tsnePerplexity.setValue(30);
    _usePreComputedTSNE.setSerializationName("CSCGDV:Use Precomputed TSNE");
    _usePreComputedTSNE.setChecked(true);
    _applyLogTransformation.setChecked(false);
    _hiddenShowncolumns.setSerializationName("CSCGDV:Hidden Shown Columns");
    _speciesExplorerInMap.setSerializationName("CSCGDV:Species Explorer In Map");
    _scatterplotReembedColorOption.setSerializationName("CSCGDV:Scatterplot Reembedding Color Option");
    _scatterplotEmbeddingPointsUMAPOption.setSerializationName("CSCGDV:Scatterplot Embedding UMAP Points Option");
    _typeofTopNGenes.setSerializationName("CSCGDV:Type of Top N Genes");
    _clusterCountSortingType.setSerializationName("CSCGDV:Cluster Count Sorting Type");
    _applyLogTransformation.setSerializationName("CSCGDV:Apply Log Transformation");
    _performGeneTableTsneAction.setChecked(false);
    _createRowMultiSelectTree.setDisabled(true);
    _selectedRowIndex.setDisabled(true);
    _selectedRowIndex.setString("");
    _scatterplotReembedColorOption.initialize({"Species","Cluster","Expression"}, "Species");
    _typeofTopNGenes.initialize({"Absolute","Negative","Positive","Mixed"}, "Positive");
    _clusterCountSortingType.initialize({ "Count","Name","Hierarchy View"}, "Count");
    _topNGenesFilter.setDefaultWidgetFlags(IntegralAction::WidgetFlag::SpinBox);

    QIcon updateIcon = Application::getIconFont("FontAwesome").getIcon("play");
    _startComputationTriggerAction.setIcon(updateIcon);
    _startComputationTriggerAction.setDefaultWidgetFlags(TriggerAction::WidgetFlag::IconText);

    QIcon exploreIcon = Application::getIconFont("FontAwesome").getIcon("wpexplorer");
    _speciesExplorerInMapTrigger.setIcon(exploreIcon);
    _speciesExplorerInMapTrigger.setDefaultWidgetFlags(TriggerAction::WidgetFlag::IconText);

    QIcon removeIcon = Application::getIconFont("FontAwesome").getIcon("backspace");
    _removeRowSelection.setIcon(removeIcon);
    _removeRowSelection.setDefaultWidgetFlags(TriggerAction::WidgetFlag::IconText);

    QIcon revertIcon = Application::getIconFont("FontAwesome").getIcon("undo");
    _revertRowSelectionChangesToInitial.setIcon(revertIcon);
    _revertRowSelectionChangesToInitial.setDefaultWidgetFlags(TriggerAction::WidgetFlag::IconText);

    _scatterplotEmbeddingPointsUMAPOption.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == PointType;
        });
    _filteringEditTreeDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == CrossSpeciesComparisonTreeType;
        });
    _referenceTreeDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == CrossSpeciesComparisonTreeType;
        });
    _mainPointsDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == PointType;
        });
    _speciesNamesDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == ClusterType;
        });
    _clusterNamesDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == ClusterType;
        });
    _embeddingDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == PointType;
        });
    const auto changetooltipTopN = [this]() -> void
        {
            _topNGenesFilter.setToolTip("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setText("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setIconText("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setToolTip("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setObjectName("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setWhatsThis("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
        };
    connect(&_topNGenesFilter, &IntegralAction::valueChanged, this, changetooltipTopN);

    const auto updatespeciesExplorerInMap = [this]() -> void
        {

            QStringList leafValues = _speciesExplorerInMap.getSelectedOptions();


            removeSelectionTableRows(&leafValues);
            enableDisableButtonsAutomatically();
        };
    connect(&_speciesExplorerInMap, &OptionsAction::selectedOptionsChanged, this, updatespeciesExplorerInMap);

    const auto updateGeneFilteringTrigger = [this]() -> void
        {
            disableActions();
            _speciesExplorerInMap.setSelectedOptions({});
            //_erroredOutFlag = false;
            updateButtonTriggered();
            enableActions();

        };
    connect(&_startComputationTriggerAction, &TriggerAction::triggered, this, updateGeneFilteringTrigger);
    const auto updateCreateRowMultiSelectTreeTrigger = [this]() -> void {
        
        if (_filteringEditTreeDataset.getCurrentDataset().isValid())
        {
            auto treeDataset = mv::data().getDataset<CrossSpeciesComparisonTree>(_filteringEditTreeDataset.getCurrentDataset().getDatasetId());

            QStringList selectedRowsStrList = _geneNamesConnection.getString().split("*%$@*@$%*");


            if (treeDataset.isValid() && selectedRowsStrList.size() > 0)
            {
                //if (speciesSelectedIndicesCounter.size() > 0)
                {
                    //QJsonObject valueStringReference = createJsonTree(speciesSelectedIndicesCounter);
                    //if (!valueStringReference.isEmpty())
                    {
                        //treeDataset->setTreeData(valueStringReference);
                        //events().notifyDatasetDataChanged(treeDataset);
                        //TODO:: add the tree to the tree dataset
/*
                        QString treeData = createJsonTreeFromNewick(QString::fromStdString(modifiedNewick), leafnames);
                        if (!treeData.isEmpty())
                        {

                            QJsonObject valueStringReference = QJsonDocument::fromJson(treeData.toUtf8()).object();
                            if (!valueStringReference.isEmpty())
                            {
                                treeDataset->setTreeData(valueStringReference);
                                events().notifyDatasetDataChanged(treeDataset);
                            }
                        }
                        */

                    }
                }


            }




        }


        };

    connect(&_createRowMultiSelectTree, &TriggerAction::triggered, this, updateCreateRowMultiSelectTreeTrigger);

    const auto updateMainPointsDataset = [this]() -> void {

        if (!_meanMapComputed)
        {
            computeGeneMeanExpressionMap();
        }
        
        
        if (_mainPointsDataset.getCurrentDataset().isValid())
        {
            
            auto fullDataset=mv::data().getDataset<Points>(_mainPointsDataset.getCurrentDataset().getDatasetId());
            auto dimensions = fullDataset->getNumDimensions();
            if (dimensions > 0) {
                _topNGenesFilter.setMinimum(0);
                _topNGenesFilter.setMaximum(dimensions);

                _topNGenesFilter.setValue(std::min(10, static_cast<int>(dimensions)));
            }
            else {
                _topNGenesFilter.setMinimum(0);
                _topNGenesFilter.setMaximum(0);
                _topNGenesFilter.setValue(0);
            }

            _topNGenesFilter.setToolTip("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setText("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setIconText("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setObjectName("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setWhatsThis("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");


        }
        else
        {
            _topNGenesFilter.setMinimum(0);
            _topNGenesFilter.setMaximum(0);
            _topNGenesFilter.setValue(0);
            _topNGenesFilter.setToolTip("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setText("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setIconText("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setObjectName("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            _topNGenesFilter.setWhatsThis("Top N genes: 0 to " + QString::number(_topNGenesFilter.getMaximum()) + " (Current: " + QString::number(_topNGenesFilter.getValue()) + ")");
            
        }
        
 };

    connect(&_mainPointsDataset, &DatasetPickerAction::currentIndexChanged, this, updateMainPointsDataset);

    const auto updateSpeciesNameDataset = [this]() -> void {
        _selectedSpeciesCellCountMap.clear();
        QStringList speciesOptions = {};
        if (_speciesNamesDataset.getCurrentDataset().isValid())
        {
            auto clusterFullDataset=mv::data().getDataset<Clusters>(_speciesNamesDataset.getCurrentDataset().getDatasetId());
            auto clusterValuesData = clusterFullDataset->getClusters();
            if (!clusterValuesData.isEmpty())
            {
                for (auto clusters : clusterValuesData)
                {
                    {
                        auto name = clusters.getName();
                        auto color = clusters.getColor();
                        _selectedSpeciesCellCountMap[name].color = color;
                        _selectedSpeciesCellCountMap[name].selectedCellsCount = 0;
                        _selectedSpeciesCellCountMap[name].nonSelectedCellsCount = 0;
                        speciesOptions.append(name);
                    }

                }

            }

        }
        _speciesExplorerInMap.setOptions(speciesOptions);
        
        if (!_meanMapComputed)
        {
            computeGeneMeanExpressionMap();
        }
        };

    connect(&_speciesNamesDataset, &DatasetPickerAction::currentIndexChanged, this, updateSpeciesNameDataset);

    const auto updateScatterplotColor = [this]() -> void {
        auto selectedColorType= _scatterplotReembedColorOption.getCurrentText();
        if (selectedColorType != "")
        {
            auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
            mv::gui::DatasetPickerAction* colorDatasetPickerAction;
            mv::gui::DatasetPickerAction* pointDatasetPickerAction;
            
            if (scatterplotViewFactory) {
                for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                    if (plugin->getGuiName() == "Scatterplot Gene Similarity View") {
                        pointDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Position"));
                        if (pointDatasetPickerAction) {


                            if(pointDatasetPickerAction->getCurrentDataset() == _selectedPointsTSNEDataset){
                            colorDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Color"));
                            if (colorDatasetPickerAction)
                            {
                                



                                auto selectedColorType = _scatterplotReembedColorOption.getCurrentText();
                                if (selectedColorType != "")
                                {
                                    if (selectedColorType == "Cluster")
                                    {
                                        if (_tsneDatasetClusterColors.isValid())
                                        {
                                            colorDatasetPickerAction->setCurrentText("");
                                            colorDatasetPickerAction->setCurrentDataset(_tsneDatasetClusterColors);
                                        }
                                    }
                                    else if (selectedColorType == "Species")
                                    {
                                        if (_tsneDatasetSpeciesColors.isValid())
                                        {
                                            colorDatasetPickerAction->setCurrentText("");
                                            colorDatasetPickerAction->setCurrentDataset(_tsneDatasetSpeciesColors);
                                        }
                                    }
                                    else if (selectedColorType == "Expression")
                                    {
                                        if (_tsneDatasetExpressionColors.isValid())
                                        {
                                            colorDatasetPickerAction->setCurrentText("");
                                            colorDatasetPickerAction->setCurrentDataset(_tsneDatasetExpressionColors);
                                        }
                                    }



                                }





                            }
                        }
                        }
                    }
                }
            }

        }
        
        };
    connect(&_scatterplotReembedColorOption, &OptionAction::currentIndexChanged, this, updateScatterplotColor);
    
    
    const auto updateStatus = [this]() -> void {
        auto string = _statusColorAction.getString();
        QString labelText = "";
        QString backgroundColor = "none";
        if (string == "C")
        {
            _startComputationTriggerAction.setDisabled(true);
        }
        else
        {
            _startComputationTriggerAction.setDisabled(false);
        }


        if (string == "C") {
            labelText = "Updated";
            backgroundColor = "#28a745"; // Green
            
        }
        else if (string == "M") {
            labelText = "Pending";
            backgroundColor = "#ffc107"; // Gold

        }
        else if (string == "E") {
            labelText ="Error Occurred! Try updating again";
            backgroundColor = "#dc3545"; // Red
        }
        else if (string == "R")
        {
            labelText = "Running updation";
            backgroundColor = "#007bff"; // Blue
        }
        else {
            labelText = "Unknown";
            backgroundColor = "#6c757d"; // Grey
        }




        // Update the _statusBarActionWidget with the new label text and background color
        _statusBarActionWidget->showMessage("Status: " + labelText);
        _statusBarActionWidget->setStyleSheet("QStatusBar{padding-left:8px;background:" + backgroundColor + ";color:white;}");


        };
    connect(&_statusColorAction, &StringAction::stringChanged, this, updateStatus);

    /*const auto updateSelectedCellClusterInfoBox = [this]() -> void {

        
        // Clear any previous message
        _selectedCellClusterInfoStatusBar->clearMessage();

        // Check if there's a previously added label and remove it
        if (_currentCellSelectionClusterInfoLabel != nullptr) {
            _selectedCellClusterInfoStatusBar->removeWidget(_currentCellSelectionClusterInfoLabel);
            delete _currentCellSelectionClusterInfoLabel; // Delete the previous label to avoid memory leaks
            _currentCellSelectionClusterInfoLabel = nullptr; // Reset the pointer to indicate there's no current label
        }

        // Create a new QLabel
        _currentCellSelectionClusterInfoLabel = new QLabel;
        auto string = _selectedCellClusterInfoBox.getString();
        QString htmlText = string; 
        _currentCellSelectionClusterInfoLabel->setText(htmlText); 
        _selectedCellClusterInfoStatusBar->addWidget(_currentCellSelectionClusterInfoLabel); 
        

        QLayoutItem* layoutItem;

        while ((layoutItem = _selectedCellClusterInfoStatusBar->takeAt(0)) != nullptr) {
            delete layoutItem->widget();
            delete layoutItem;
        }

        for (cluster : clusters) {
            auto clusterLabel = new QLabel(parent, clusterName);
            clusterLabel->setStyleSheet("");
            _selectedCellClusterInfoStatusBar->addWidget(clusterLabel);
        }



        };*/

    const auto updateEmbeddingDataset = [this]() -> void {


        };
    connect(&_embeddingDataset, &DatasetPickerAction::currentIndexChanged, this, updateEmbeddingDataset);


    const auto updateTypeOfTopNGenesFilter = [this]() -> void {
        _statusColorAction.setString("M");

        };
    connect(&_typeofTopNGenes, &OptionAction::currentIndexChanged, this, updateTypeOfTopNGenesFilter);

    const auto updateApplyLogTransformation = [this]() -> void {
        _statusColorAction.setString("M");

        };
    connect(&_applyLogTransformation, &ToggleAction::toggled, this, updateApplyLogTransformation);

    const auto updateClusterCountSortingType = [this]() -> void {
        _statusColorAction.setString("M");

        };
    connect(&_clusterCountSortingType, &OptionAction::currentIndexChanged, this, updateClusterCountSortingType);

    const auto updateTopGenesSlider = [this]() -> void {
        _statusColorAction.setString("M");

        };
    connect(&_topNGenesFilter, &IntegralAction::valueChanged, this, updateTopGenesSlider);

    _statusColorAction.setString("M");
}



void SettingsAction::updateButtonTriggered()
{
    try {
       // _startComputationTriggerAction.setDisabled(true);
        startCodeTimer("UpdateGeneFilteringTrigger");
        //startCodeTimer("Part1");


        auto pointsDataset = _mainPointsDataset.getCurrentDataset();
        auto embeddingDataset = _embeddingDataset.getCurrentDataset();
        auto speciesDataset = _speciesNamesDataset.getCurrentDataset();
        auto clusterDataset = _clusterNamesDataset.getCurrentDataset();
        auto referenceTreeDataset = _referenceTreeDataset.getCurrentDataset();
        _selectedSpeciesVals.setString("");
        _geneNamesConnection.setString("");
        bool isValid = false;

        QString referenceTreedatasetId = "";
        //stopCodeTimer("Part1");
        //startCodeTimer("Part2");
        if (!pointsDataset.isValid() || !embeddingDataset.isValid() || !speciesDataset.isValid() || !clusterDataset.isValid() || !referenceTreeDataset.isValid())
        {
            qDebug() << "No datasets selected";
            //_startComputationTriggerAction.setDisabled(false);
            return;
        }
        if (pointsDataset->getSelectionIndices().size() < 1)
        {
            qDebug() << "No points selected";
            //_startComputationTriggerAction.setDisabled(false);
            return;
        }
        if (_selectedPointsTSNEDataset.isValid())
        {
            _selectedPointsTSNEDataset->setSelectionIndices({});
        }
        //stopCodeTimer("Part2");
        //startCodeTimer("Part3");
        _clusterNameToGeneNameToExpressionValue.clear();
        referenceTreedatasetId = referenceTreeDataset->getId();
        isValid = speciesDataset->getParent() == pointsDataset && clusterDataset->getParent() == pointsDataset && embeddingDataset->getParent() == pointsDataset;
        if (!isValid)
        {
            qDebug() << "Datasets are not valid";
            //_startComputationTriggerAction.setDisabled(false);
            return;
        }
        _selectedIndicesFromStorage.clear();
        _selectedIndicesFromStorage = pointsDataset->getSelectionIndices();

        auto embeddingDatasetRaw = mv::data().getDataset<Points>(embeddingDataset->getId());
        auto pointsDatasetRaw = mv::data().getDataset<Points>(pointsDataset->getId());
        auto pointsDatasetallColumnNameList = pointsDatasetRaw->getDimensionNames();
        auto embeddingDatasetallColumnNameList = embeddingDatasetRaw->getDimensionNames();
        //stopCodeTimer("Part3");
        //startCodeTimer("Part4");
        std::vector<int> embeddingDatasetColumnIndices(embeddingDatasetallColumnNameList.size());
        std::iota(embeddingDatasetColumnIndices.begin(), embeddingDatasetColumnIndices.end(), 0);

        std::vector<int> pointsDatasetallColumnIndices(pointsDatasetallColumnNameList.size());
        std::iota(pointsDatasetallColumnIndices.begin(), pointsDatasetallColumnIndices.end(), 0);
        //stopCodeTimer("Part4");
        {

            if (_selectedIndicesFromStorage.size() > 0 && embeddingDatasetColumnIndices.size() > 0)
            {
                //startCodeTimer("Part5");
                auto speciesDatasetRaw = mv::data().getDataset<Clusters>(speciesDataset->getId());
                auto clusterDatasetRaw = mv::data().getDataset<Clusters>(clusterDataset->getId());
                auto clusterDatasetName = clusterDatasetRaw->getGuiName();
                auto clustersValuesAll = clusterDatasetRaw->getClusters();
                auto speciesValuesAll = speciesDatasetRaw->getClusters();

                std::map<QString, std::pair<QColor, std::vector<int>>> selectedClustersMap;
                std::map<QString, std::pair<QColor, std::vector<int>>> selectedSpeciesMap;
                //stopCodeTimer("Part5");
                if (!speciesValuesAll.empty() && !clustersValuesAll.empty())
                {
                    //startCodeTimer("Part6.1");
                    if (!_selectedPointsDataset.isValid())
                    {
                        _selectedPointsDataset = mv::data().createDataset("Points", "SelectedPointsDataset");
                        _selectedPointsDataset->setGroupIndex(10);
                        mv::events().notifyDatasetAdded(_selectedPointsDataset);

                    }
                    if (!_tsneDatasetExpressionColors.isValid())
                    {
                        _tsneDatasetExpressionColors = mv::data().createDataset("Points", "TSNEDatasetExpressionColors", _selectedPointsDataset);
                        _tsneDatasetExpressionColors->setGroupIndex(10);
                        mv::events().notifyDatasetAdded(_tsneDatasetExpressionColors);

                    }
                    if (!_filteredUMAPDatasetPoints.isValid())
                    {
                        _filteredUMAPDatasetPoints = mv::data().createDataset("Points", "Filtered UMAP Dataset Points");
                        _filteredUMAPDatasetPoints->setGroupIndex(10);
                        mv::events().notifyDatasetAdded(_filteredUMAPDatasetPoints);
                        if (!_filteredUMAPDatasetColors.isValid())
                        {
                            //need to delete

                        }
                        _filteredUMAPDatasetColors = mv::data().createDataset("Points", "Filtered UMAP Dataset Colors", _filteredUMAPDatasetPoints);
                        _filteredUMAPDatasetColors->setGroupIndex(10);
                        mv::events().notifyDatasetAdded(_filteredUMAPDatasetColors);

                    }
                    if (!_selectedPointsEmbeddingDataset.isValid())
                    {
                        _selectedPointsEmbeddingDataset = mv::data().createDataset("Points", "TSNEDataset", _selectedPointsDataset);
                        _selectedPointsEmbeddingDataset->setGroupIndex(10);
                        mv::events().notifyDatasetAdded(_selectedPointsEmbeddingDataset);

                    }
                    if (!_tsneDatasetSpeciesColors.isValid())
                    {
                        _tsneDatasetSpeciesColors = mv::data().createDataset("Cluster", "TSNEDatasetSpeciesColors", _selectedPointsDataset);
                        _tsneDatasetSpeciesColors->setGroupIndex(10);
                        mv::events().notifyDatasetAdded(_tsneDatasetSpeciesColors);
                    }

                    if (!_tsneDatasetClusterColors.isValid())
                    {
                        _tsneDatasetClusterColors = mv::data().createDataset("Cluster", "TSNEDatasetClusterColors", _selectedPointsDataset);
                        _tsneDatasetClusterColors->setGroupIndex(10);
                        mv::events().notifyDatasetAdded(_tsneDatasetClusterColors);
                    }
                    //stopCodeTimer("Part6.1");
                    if (_selectedPointsDataset.isValid() && _selectedPointsEmbeddingDataset.isValid() && _tsneDatasetSpeciesColors.isValid() && _tsneDatasetClusterColors.isValid())
                    {
                        //startCodeTimer("Part6.2");
                        _tsneDatasetSpeciesColors->getClusters() = QVector<Cluster>();
                        events().notifyDatasetDataChanged(_tsneDatasetSpeciesColors);
                        _tsneDatasetClusterColors->getClusters() = QVector<Cluster>();
                        events().notifyDatasetDataChanged(_tsneDatasetClusterColors);
                        //stopCodeTimer("Part6.2");
                        //startCodeTimer("Part7");
                        //startCodeTimer("Part7.1");
                        int selectedIndicesFromStorageSize = _selectedIndicesFromStorage.size();
                        int pointsDatasetColumnsSize = pointsDatasetallColumnIndices.size();
                        int embeddingDatasetColumnsSize = embeddingDatasetColumnIndices.size();
                        QString datasetIdEmb = _selectedPointsDataset->getId();
                        QString datasetId = _selectedPointsEmbeddingDataset->getId();
                        int dimofDatasetExp = 1;
                        std::vector<QString> dimensionNamesExp = { "Expression" };
                        QString datasetIdExp = _tsneDatasetExpressionColors->getId();
                        //stopCodeTimer("Part7.1");
                        //startCodeTimer("Part7.2");

                        // Define result containers outside the lambda functions to ensure they are accessible later
                        std::vector<float> resultContainerForSelectedPoints(selectedIndicesFromStorageSize * pointsDatasetColumnsSize);
                        std::vector<float> resultContainerForSelectedEmbeddingPoints(selectedIndicesFromStorageSize * embeddingDatasetColumnsSize);
                        std::vector<float> resultContainerColorPoints(selectedIndicesFromStorageSize, -1.0f);

                        //first thread start
                        auto future1 = std::async(std::launch::async, [&]() {
                            pointsDatasetRaw->populateDataForDimensions(resultContainerForSelectedPoints, pointsDatasetallColumnIndices, _selectedIndicesFromStorage);
                            });

                        //second thread start
                        auto future2 = std::async(std::launch::async, [&]() {
                            embeddingDatasetRaw->populateDataForDimensions(resultContainerForSelectedEmbeddingPoints, embeddingDatasetColumnIndices, _selectedIndicesFromStorage);
                            });


                        // Wait for all futures to complete before proceeding
                        future1.wait();
                        future2.wait();


                        //startCodeTimer("Part7.2.1");
                        //needs to wait for future1 finish only
                        populatePointData(datasetIdEmb, resultContainerForSelectedPoints, selectedIndicesFromStorageSize, pointsDatasetColumnsSize, pointsDatasetallColumnNameList);
                        //stopCodeTimer("Part7.2.1");

                        //startCodeTimer("Part7.2.2");
                        //needs to wait for future2 finish only
                        populatePointData(datasetId, resultContainerForSelectedEmbeddingPoints, selectedIndicesFromStorageSize, embeddingDatasetColumnsSize, embeddingDatasetallColumnNameList);
                        //stopCodeTimer("Part7.2.2");

                        //startCodeTimer("Part7.2.3");
                        //needs to wait for future3 finish only
                        populatePointData(datasetIdExp, resultContainerColorPoints, selectedIndicesFromStorageSize, dimofDatasetExp, dimensionNamesExp);
                        //stopCodeTimer("Part7.2.3");

                        //stopCodeTimer("Part7.2");


                        //stopCodeTimer("Part7");
                        //startCodeTimer("Part8");
                        if (_selectedPointsTSNEDataset.isValid())
                        {
                            auto runningAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDataset->findChildByPath("TSNE/TsneComputationAction/Running"));

                            if (runningAction)
                            {

                                if (runningAction->isChecked())
                                {
                                    auto stopAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDataset->findChildByPath("TSNE/TsneComputationAction/Stop"));
                                    if (stopAction)
                                    {
                                        stopAction->trigger();
                                        std::this_thread::sleep_for(std::chrono::seconds(5));
                                    }
                                }

                            }
                        }

                        if (_selectedPointsTSNEDataset.isValid())
                        {
                            auto datasetIDLowRem = _selectedPointsTSNEDataset.getDatasetId();
                            mv::events().notifyDatasetAboutToBeRemoved(_selectedPointsTSNEDataset);
                            mv::data().removeDataset(_selectedPointsTSNEDataset);
                            mv::events().notifyDatasetRemoved(datasetIDLowRem, PointType);
                        }
                        //stopCodeTimer("Part8");
                        //startCodeTimer("Part9");
                        mv::plugin::AnalysisPlugin* analysisPlugin;
                        bool usePreTSNE = _usePreComputedTSNE.isChecked();

                        auto scatterplotModificationsLowDimUMAP = [this]() {
                            if (_selectedPointsTSNEDataset.isValid()) {
                                auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
                                mv::gui::DatasetPickerAction* colorDatasetPickerAction;
                                mv::gui::DatasetPickerAction* pointDatasetPickerAction;
                                if (scatterplotViewFactory) {
                                    for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                                        if (plugin->getGuiName() == "Scatterplot Gene Similarity View") {
                                            pointDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Position"));
                                            if (pointDatasetPickerAction) {
                                                pointDatasetPickerAction->setCurrentText("");

                                                pointDatasetPickerAction->setCurrentDataset(_selectedPointsTSNEDataset);

                                                colorDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Color"));
                                                if (colorDatasetPickerAction)
                                                {
                                                    colorDatasetPickerAction->setCurrentText("");



                                                    auto selectedColorType = _scatterplotReembedColorOption.getCurrentText();
                                                    if (selectedColorType != "")
                                                    {
                                                        if (selectedColorType == "Cluster")
                                                        {
                                                            if (_tsneDatasetClusterColors.isValid())
                                                            {
                                                                colorDatasetPickerAction->setCurrentDataset(_tsneDatasetClusterColors);
                                                            }
                                                        }
                                                        else if (selectedColorType == "Species")
                                                        {
                                                            if (_tsneDatasetSpeciesColors.isValid())
                                                            {
                                                                colorDatasetPickerAction->setCurrentDataset(_tsneDatasetSpeciesColors);
                                                            }
                                                        }
                                                        else if (selectedColorType == "Expression")
                                                        {
                                                            if (_tsneDatasetExpressionColors.isValid())
                                                            {
                                                                colorDatasetPickerAction->setCurrentDataset(_tsneDatasetExpressionColors);
                                                            }
                                                        }



                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            };
                        //stopCodeTimer("Part9");
                        if (!usePreTSNE)
                        {
                            //startCodeTimer("Part10");
                            analysisPlugin = mv::plugins().requestPlugin<AnalysisPlugin>("tSNE Analysis", { _selectedPointsEmbeddingDataset });
                            if (!analysisPlugin) {
                                qDebug() << "Could not find create TSNE Analysis";
                                return;
                            }
                            _selectedPointsTSNEDataset = analysisPlugin->getOutputDataset();
                            if (_selectedPointsTSNEDataset.isValid())
                            {

                                int perplexity = std::min(static_cast<int>(_selectedIndicesFromStorage.size()), _tsnePerplexity.getValue());
                                if (perplexity < 5)
                                {
                                    qDebug() << "Perplexity is less than 5";
                                    //_startComputationTriggerAction.setDisabled(false);
                                    return;
                                }
                                if (perplexity != _tsnePerplexity.getValue())
                                {
                                    _tsnePerplexity.setValue(perplexity);
                                }

                                auto perplexityAction = dynamic_cast<IntegralAction*>(_selectedPointsTSNEDataset->findChildByPath("TSNE/Perplexity"));
                                if (perplexityAction)
                                {
                                    qDebug() << "Perplexity: Found";
                                    perplexityAction->setValue(perplexity);
                                }
                                else
                                {
                                    qDebug() << "Perplexity: Not Found";
                                }

                                scatterplotModificationsLowDimUMAP();

                                auto startAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDataset->findChildByPath("TSNE/TsneComputationAction/Start"));
                                if (startAction) {

                                    startAction->trigger();

                                    analysisPlugin->getOutputDataset()->setSelectionIndices({});
                                }

                            }
                            //stopCodeTimer("Part10");
                        }
                        else
                        {
                            //startCodeTimer("Part11");
                            auto umapDataset = _scatterplotEmbeddingPointsUMAPOption.getCurrentDataset();

                            if (umapDataset.isValid())
                            {


                                _selectedPointsTSNEDataset = mv::data().createDerivedDataset<Points>("SelectedPointsTSNEDataset", _selectedPointsEmbeddingDataset, _selectedPointsEmbeddingDataset);

                                mv::events().notifyDatasetAdded(_selectedPointsTSNEDataset);

                                auto umapDatasetRaw = mv::data().getDataset<Points>(umapDataset->getId());
                                auto dimNames = umapDatasetRaw->getDimensionNames();
                                int preComputedEmbeddingColumnsSize = umapDatasetRaw->getNumDimensions();
                                std::vector<float> resultContainerPreComputedUMAP(selectedIndicesFromStorageSize * preComputedEmbeddingColumnsSize);
                                std::vector<int> preComputedEmbeddingColumnIndices(preComputedEmbeddingColumnsSize);

                                std::iota(preComputedEmbeddingColumnIndices.begin(), preComputedEmbeddingColumnIndices.end(), 0);

                                umapDatasetRaw->populateDataForDimensions(resultContainerPreComputedUMAP, preComputedEmbeddingColumnIndices, _selectedIndicesFromStorage);

                                QString datasetId = _selectedPointsTSNEDataset->getId();
                                populatePointData(datasetId, resultContainerPreComputedUMAP, selectedIndicesFromStorageSize, preComputedEmbeddingColumnsSize, dimNames);

                                if (_selectedPointsTSNEDataset.isValid())
                                {
                                    scatterplotModificationsLowDimUMAP();
                                }
                            }
                            else
                            {
                                qDebug() << "UMAP Dataset not valid";
                            }


                            //stopCodeTimer("Part11");



                        }



                    }
                    else
                    {
                        qDebug() << "Datasets are not valid";
                    }
                    //startCodeTimer("Part12");
                    //startCodeTimer("Part12.1");
                    QFuture<void> futureClusterCVals = QtConcurrent::run([&]() {
                        QMutex mutex;
                        for (auto& clusters : clustersValuesAll) {
                            auto clusterIndices = clusters.getIndices();
                            auto clusterName = clusters.getName();
                            auto clusterColor = clusters.getColor();
                            std::vector<int> filteredIndices;

                            QtConcurrent::blockingMap(clusterIndices, [&](int index) {
                                int indexVal = findIndex(_selectedIndicesFromStorage, index);
                                if (indexVal != -1) {
                                    QMutexLocker locker(&mutex);
                                    filteredIndices.push_back(indexVal);
                                }
                                });

                            {
                                QMutexLocker locker(&mutex);
                                selectedClustersMap[clusterName] = { clusterColor, filteredIndices };
                            }
                        }
                        });
                    QFuture<void> futureSpeciesCVals = QtConcurrent::run([&]() {
                        QMutex mutex;
                        for (auto& clusters : speciesValuesAll) {
                            auto clusterIndices = clusters.getIndices();
                            auto clusterName = clusters.getName();
                            auto clusterColor = clusters.getColor();
                            std::vector<int> filteredIndices;

                            QtConcurrent::blockingMap(clusterIndices, [&](int index) {
                                int indexVal = findIndex(_selectedIndicesFromStorage, index);
                                if (indexVal != -1) {
                                    QMutexLocker locker(&mutex);
                                    filteredIndices.push_back(indexVal);
                                }
                                });

                            {
                                QMutexLocker locker(&mutex);
                                selectedSpeciesMap[clusterName] = { clusterColor, filteredIndices };
                            }
                        }
                        });
                    futureClusterCVals.waitForFinished(); // Wait for the concurrent task to complete
                    futureSpeciesCVals.waitForFinished();
                    //stopCodeTimer("Part12.1");

                    //startCodeTimer("Part12.2");


                    std::sort(_selectedIndicesFromStorage.begin(), _selectedIndicesFromStorage.end());

                    QMutex clusterGeneMeanExpressionMapMutex;
                    QMutex clusterNameToGeneNameToExpressionValueMutex;

                    QFutureSynchronizer<void> synchronizer; // Use QFutureSynchronizer to wait for all tasks

                    for (auto& species : speciesValuesAll) {
                        auto future = QtConcurrent::run([&, species]() {
                            auto speciesIndices = species.getIndices();
                            auto speciesName = species.getName();

                            std::vector<int> commonSelectedIndices;

                            std::sort(speciesIndices.begin(), speciesIndices.end());
                            std::set_intersection(_selectedIndicesFromStorage.begin(), _selectedIndicesFromStorage.end(), speciesIndices.begin(), speciesIndices.end(), std::back_inserter(commonSelectedIndices));

                            std::unordered_map<QString, Stats> localClusterNameToGeneNameToExpressionValue;

                            for (int i = 0; i < pointsDatasetallColumnNameList.size(); i++) {
                                const auto& geneName = pointsDatasetallColumnNameList[i];
                                std::vector<int> geneIndex = { i };

                                // Access shared data in a thread-safe manner
                                QMutexLocker locker(&clusterGeneMeanExpressionMapMutex);
                                const auto& nonSelectionDetails = _clusterGeneMeanExpressionMap[speciesName][geneName];
                                locker.unlock();

                                int allCellCounts = nonSelectionDetails.first;
                                float allCellMean = nonSelectionDetails.second;

                                float nonSelectedMean = 0.0;
                                int nonSelectedCells = 0;

                                StatisticsSingle calculateStatisticsShort = { 0.0f, 0 };

                                if (!commonSelectedIndices.empty()) {
                                    std::vector<float> resultContainerShort(commonSelectedIndices.size());
                                    pointsDatasetRaw->populateDataForDimensions(resultContainerShort, geneIndex, commonSelectedIndices);
                                    calculateStatisticsShort = calculateStatistics(resultContainerShort);

                                    float allCellTotal = allCellMean * allCellCounts;

                                    nonSelectedCells = allCellCounts - calculateStatisticsShort.countVal;

                                    // Check to prevent division by zero
                                    if (nonSelectedCells > 0) {
                                        nonSelectedMean = (allCellTotal - calculateStatisticsShort.meanVal * calculateStatisticsShort.countVal) / nonSelectedCells;
                                    }
                                    else {
                                        nonSelectedMean = 0.0f;
                                    }
                                }
                                else {
                                    nonSelectedMean = allCellMean;
                                    nonSelectedCells = allCellCounts;
                                }

                                StatisticsSingle calculateStatisticsNot = { nonSelectedMean, nonSelectedCells };

                                localClusterNameToGeneNameToExpressionValue[geneName] = combineStatisticsSingle(calculateStatisticsShort, calculateStatisticsNot);
                            }

                            // Merge results in a thread-safe manner
                            QMutexLocker locker(&clusterNameToGeneNameToExpressionValueMutex);
                            for (const auto& pair : localClusterNameToGeneNameToExpressionValue) {
                                _clusterNameToGeneNameToExpressionValue[speciesName][pair.first] = pair.second;
                            }
                            });
                        synchronizer.addFuture(future); // Add each future to the synchronizer
                    }

                    synchronizer.waitForFinished();



                    //stopCodeTimer("Part12.2");

                    auto clusterColorDatasetId = _tsneDatasetClusterColors->getId();
                    auto speciesColorDatasetId = _tsneDatasetSpeciesColors->getId();
                    //startCodeTimer("Part12.3");
                    populateClusterData(speciesColorDatasetId, selectedSpeciesMap);
                    //stopCodeTimer("Part12.3");
                    //startCodeTimer("Part12.4");
                    populateClusterData(clusterColorDatasetId, selectedClustersMap);
                    //stopCodeTimer("Part12.4");
                    //stopCodeTimer("Part12");
                    if (_tsneDatasetClusterColors.isValid())
                    {

                        auto clusterValues = _tsneDatasetClusterColors->getClusters();
                        if (!clusterValues.empty())
                        {
                            //startCodeTimer("Part13");

                            QLayoutItem* layoutItem;
                            while ((layoutItem = _selectedCellClusterInfoStatusBar->takeAt(0)) != nullptr) {
                                delete layoutItem->widget();
                                delete layoutItem;
                            }

                            // Create a description label
                            auto descriptionLabel = new QLabel("Selected Cell Counts per " + clusterDatasetName + " :");
                            // Optionally, set a stylesheet for the description label for styling
                            descriptionLabel->setStyleSheet("QLabel { font-weight: bold; padding: 2px; }");
                            // Add the description label to the layout
                            _selectedCellClusterInfoStatusBar->addWidget(descriptionLabel);


                            std::vector<ClusterOrderContainer> orderedClustersSet;

                            for (const auto& cluster : clusterValues) {
                                ClusterOrderContainer temp{
                                    cluster.getIndices().size(),
                                    cluster.getColor(),
                                    cluster.getName()
                                };
                                orderedClustersSet.push_back(std::move(temp));
                            }

                            const auto& currentText = _clusterCountSortingType.getCurrentText();
                            if (currentText == "Name") {
                                std::sort(orderedClustersSet.begin(), orderedClustersSet.end(), sortByName);
                            }
                            else if (currentText == "Hierarchy View" && !_customOrderClustersFromHierarchy.empty()) {
                                if (_customOrderClustersFromHierarchyMap.empty()) {
                                    _customOrderClustersFromHierarchyMap = prepareCustomSortMap(_customOrderClustersFromHierarchy);
                                }
                                std::sort(orderedClustersSet.begin(), orderedClustersSet.end(), [&](const ClusterOrderContainer& a, const ClusterOrderContainer& b) {
                                    return sortByCustomList(a, b, _customOrderClustersFromHierarchyMap);
                                    });
                            }
                            else {
                                std::sort(orderedClustersSet.begin(), orderedClustersSet.end(), sortByCount);
                                if (currentText != "Count") {
                                    _clusterCountSortingType.setCurrentText("Count");
                                }
                            }

                            for (const auto& clustersFromSet : orderedClustersSet)
                            {
                                auto clusterLabel = new QLabel(QString("%1: %2").arg(clustersFromSet.name).arg(clustersFromSet.count));
                                QColor textColor = clustersFromSet.color.lightness() > 127 ? Qt::black : Qt::white;
                                clusterLabel->setStyleSheet(QString("QLabel { color: %1; background-color: %2; padding: 2px; border: 0.5px solid %3; }")
                                    .arg(textColor.name()).arg(clustersFromSet.color.name(QColor::HexArgb)).arg(textColor.name()));
                                _selectedCellClusterInfoStatusBar->addWidget(clusterLabel);
                            }



                            /*

                            for (auto cluster : clusterValues) {
                                auto clusterName = cluster.getName();
                                auto clusterIndicesSize = cluster.getIndices().size();
                                auto clusterColor = cluster.getColor(); // Assuming getColor() returns a QColor

                                // Calculate luminance
                                qreal luminance = 0.299 * clusterColor.redF() + 0.587 * clusterColor.greenF() + 0.114 * clusterColor.blueF();

                                // Choose text color based on luminance
                                QString textColor = (luminance > 0.5) ? "black" : "white";

                                // Convert QColor to hex string for stylesheet
                                QString backgroundColor = clusterColor.name(QColor::HexArgb);

                                auto clusterLabel = new QLabel(QString("%1: %2").arg(clusterName).arg(clusterIndicesSize));
                                // Add text color and background color to clusterLabel with padding and border for better styling
                                clusterLabel->setStyleSheet(QString("QLabel { color: %1; background-color: %2; padding: 2px; border: 0.5px solid %3; }")
                                    .arg(textColor).arg(backgroundColor).arg(textColor));
                                _selectedCellClusterInfoStatusBar->addWidget(clusterLabel);
                            }
                            */

                        }

                    }

                    //the next line should only execute if all above are finished


                    //startCodeTimer("Part14");
                    QVariant geneListTable = findTopNGenesPerCluster(_clusterNameToGeneNameToExpressionValue, _topNGenesFilter.getValue(), referenceTreedatasetId, 1.0);
                    //stopCodeTimer("Part14");

                    setModifiedTriggeredData(geneListTable);


                }

                else
                {
                    qDebug() << "Species or Clusters are empty";
                }


            }

            else
            {
                qDebug() << "No points selected or no dimensions present";
            }

            _removeRowSelection.trigger();
            _removeRowSelection.setEnabled(false);
            //enableDisableButtonsAutomatically();

        }
        stopCodeTimer("UpdateGeneFilteringTrigger");
        //_startComputationTriggerAction.setDisabled(false);
    }
    catch (const std::exception& e) {
        qDebug() << "An exception occurred in coputation: " << e.what();
    }
    catch (...) {
        qDebug() << "An unknown exception occurred in coputation";
    }
}

void SettingsAction::setModifiedTriggeredData(QVariant geneListTable)
{
    if (!geneListTable.isNull())
    {
        //startCodeTimer("Part15");
        //_filteredGeneNamesVariant.setVariant(geneListTable);
        _listModel.setVariant(geneListTable);
        //stopCodeTimer("Part15");

    }
    else
    {
        qDebug() << "QVariant empty";
    }
}

void SettingsAction::computeGeneMeanExpressionMap()
{
    
    
    _clusterGeneMeanExpressionMap.clear(); 
    if (_speciesNamesDataset.getCurrentDataset().isValid() && _mainPointsDataset.getCurrentDataset().isValid()) {
        auto speciesClusterDatasetFull = mv::data().getDataset<Clusters>(_speciesNamesDataset.getCurrentDataset().getDatasetId());
        auto mainPointDatasetFull = mv::data().getDataset<Points>(_mainPointsDataset.getCurrentDataset().getDatasetId());
        if (speciesClusterDatasetFull.isValid() && mainPointDatasetFull.isValid()) {
            auto speciesclusters = speciesClusterDatasetFull->getClusters();
            auto mainPointDimensionNames = mainPointDatasetFull->getDimensionNames();
            for (auto species : speciesclusters) {
                auto speciesIndices = species.getIndices();
                auto speciesName = species.getName();
                for (int i = 0; i < mainPointDimensionNames.size(); i++) {
                    auto& geneName = mainPointDimensionNames[i];
                    auto geneIndex = { i };
                    std::vector<float> resultContainerFull(speciesIndices.size());
                    mainPointDatasetFull->populateDataForDimensions(resultContainerFull, geneIndex, speciesIndices);
                    float fullMean = calculateMean(resultContainerFull);
                    _clusterGeneMeanExpressionMap[speciesName][geneName] = std::make_pair(speciesIndices.size(), fullMean);
                }

            }
            _meanMapComputed = true;

        }
    }

}

QVariant SettingsAction::findTopNGenesPerCluster(const std::map<QString, std::map<QString, Stats>>& map, int n, QString datasetId, float treeSimilarityScore) {
    
    if (map.empty() || n <= 0) {
        return QVariant();
    }
    startCodeTimer("findTopNGenesPerCluster");
    enum class SelectionOption {
        AbsoluteTopN,
        PositiveTopN,
        NegativeTopN,
        MixedTopN
    };
    auto optionValue = _typeofTopNGenes.getCurrentText();
    SelectionOption option = SelectionOption::AbsoluteTopN;
    if (optionValue == "Positive") {
        option = SelectionOption::PositiveTopN;
    }
    else if (optionValue == "Negative") {
        option = SelectionOption::NegativeTopN;
    }
    else if (optionValue == "Mixed") {
        option = SelectionOption::MixedTopN;
    }

    QSet<QString> geneList;
    std::map<QString, std::vector<QString>> geneAppearanceCounter;
    std::map<QString, std::vector<std::pair<QString, int>>> rankingMap;

    for (const auto& outerPair : map) {
        auto speciesName = outerPair.first;
        std::vector<std::pair<QString, float>> geneExpressionVec;
        geneExpressionVec.reserve(outerPair.second.size());
        for (const auto& innerPair : outerPair.second) {
            auto geneName = innerPair.first;
            auto differenceMeanValue = innerPair.second.meanSelected - innerPair.second.meanNonSelected;
            geneExpressionVec.push_back(std::make_pair(geneName, differenceMeanValue));
        }

        // Sort the geneExpressionVec based on the mean value from highest to lowest
        std::sort(geneExpressionVec.begin(), geneExpressionVec.end(), [](const auto& a, const auto& b) {
            return a.second > b.second;
            });

        switch (option) {
        case SelectionOption::AbsoluteTopN: {
            std::sort(geneExpressionVec.begin(), geneExpressionVec.end(), [](const auto& a, const auto& b) {
                return std::abs(a.second) > std::abs(b.second);
                });
            for (int i = 0; i < geneExpressionVec.size(); ++i) {
                if (i < n) {
                    geneList.insert(geneExpressionVec[i].first);
                    //check if the selected counter in the mpa is greater than 0
                    if (outerPair.second.find(geneExpressionVec[i].first)->second.meanSelected > 0) {
                        geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                    }

                    //geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                }

                rankingMap[geneExpressionVec[i].first].emplace_back(speciesName, i + 1);

            }


            break;
        }
        case SelectionOption::PositiveTopN: {
            for (int i = 0; i < geneExpressionVec.size(); ++i) {
                if (i < n) {
                    geneList.insert(geneExpressionVec[i].first);
                    if (outerPair.second.find(geneExpressionVec[i].first)->second.meanSelected > 0) {
                        geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                    }

                    //geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                }
                rankingMap[geneExpressionVec[i].first].emplace_back(speciesName, i + 1); // Adding rank, incremented by 1
            }
            break;
        }
        case SelectionOption::NegativeTopN: {
            std::reverse(geneExpressionVec.begin(), geneExpressionVec.end()); // Reverse to get the lowest values
            for (int i = 0; i < geneExpressionVec.size(); ++i) {
                if (i < n) {
                    geneList.insert(geneExpressionVec[i].first);
                    if (outerPair.second.find(geneExpressionVec[i].first)->second.meanSelected > 0) {
                        geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                    }

                    //geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                }
                rankingMap[geneExpressionVec[i].first].emplace_back(speciesName, geneExpressionVec.size() - i); // Corrected rank calculation
            }
            break;
        }
        case SelectionOption::MixedTopN: {
            int halfN = n / 2;
            // Process top halfN genes
            for (int i = 0; i < halfN; ++i) {
                geneList.insert(geneExpressionVec[i].first);
                if (outerPair.second.find(geneExpressionVec[i].first)->second.meanSelected > 0) {
                    geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                }

                //geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                rankingMap[geneExpressionVec[i].first].emplace_back(speciesName, i + 1);
            }
            // Process bottom halfN genes, if n is odd, include the middle element
            int startIdx = std::max(static_cast<int>(geneExpressionVec.size()) - halfN, halfN + (n % 2));
            for (int i = startIdx; i < geneExpressionVec.size(); ++i) {
                geneList.insert(geneExpressionVec[i].first);
                if (outerPair.second.find(geneExpressionVec[i].first)->second.meanSelected > 0) {
                    geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                }

                //geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                rankingMap[geneExpressionVec[i].first].emplace_back(speciesName, i + 1);
            }
            break;
        }

        }
    }



    stopCodeTimer("findTopNGenesPerCluster");
    QVariant returnedmodel= createModelFromData(geneList, map, datasetId, treeSimilarityScore, geneAppearanceCounter, rankingMap, n);
    return returnedmodel;
}

QVariant SettingsAction::createModelFromData(const QSet<QString>& returnGeneList, const std::map<QString, std::map<QString, Stats>>& map, const QString& treeDatasetId, const float& treeSimilarityScore, const std::map<QString, std::vector<QString>>& geneCounter, const std::map<QString, std::vector<std::pair<QString, int>>>& rankingMap, const int& n) {

    if (returnGeneList.isEmpty() || map.empty()) {
        return QVariant();
    }
    startCodeTimer("createModelFromData");
    QStandardItemModel* model = new QStandardItemModel();
    _initColumnNames = { "ID", "Species \nAppearance", "Gene Appearance Species Names", "Statistics" };
    model->setHorizontalHeaderLabels(_initColumnNames);

    QStringList headers = _initColumnNames;
    _hiddenShowncolumns.setOptions(headers);
    _hiddenShowncolumns.setSelectedOptions({ headers[0], headers[1] });

    for (const auto& gene : returnGeneList) {
        QList<QStandardItem*> row;
        std::vector<float> numbers;
        numbers.reserve(map.size()); // Reserve capacity based on map size
        std::map<QString, Stats> statisticsValuesForSpeciesMap;

        for (const auto& [speciesName, innerMap] : map) {
            auto it = innerMap.find(gene);
            if (it != innerMap.end()) {
                float value = it->second.meanSelected;
                numbers.push_back(value);
                statisticsValuesForSpeciesMap[speciesName] = it->second;
            }
        }

        row.push_back(new QStandardItem(gene)); // ID

        std::map<QString, int> rankcounter;
        if (auto rankit = rankingMap.find(gene); rankit != rankingMap.end()) {
            // Assuming rankit->second is of type std::vector<std::pair<QString,int>>
            for (const auto& pair : rankit->second) {
                rankcounter[pair.first] = pair.second;
            }
        }

        QString speciesGeneAppearancesComb;
        int count = 0;
        if (auto it = geneCounter.find(gene); it != geneCounter.end()) {
            const auto& speciesDetails = it->second;
            count = speciesDetails.size();
            QStringList speciesNames;
            for (const auto& speciesDetail : speciesDetails) {
                speciesNames << speciesDetail;
            }
            speciesGeneAppearancesComb = speciesNames.join(";");
        }

        row.push_back(new QStandardItem(QString::number(count))); // Gene Appearances
        row.push_back(new QStandardItem(speciesGeneAppearancesComb)); // Gene Appearance Species Names

        QString formattedStatistics;
        for (const auto& [species, stats] : statisticsValuesForSpeciesMap) {
            formattedStatistics += QString("Species: %1, Rank: %2, MeanSelected: %3, CountSelected: %4, MeanNotSelected: %5, CountNotSelected: %6;\n")//, MeanAll: %7, CountAll: %8
                .arg(species)
                .arg(rankcounter[species])
                .arg(stats.meanSelected, 0, 'f', 2)
                .arg(stats.countSelected)
                .arg(stats.meanNonSelected, 0, 'f', 2)
                .arg(stats.countNonSelected)
                //.arg(stats.meanAll, 0, 'f', 2)
                //.arg(stats.countAll)
                ;
        }
        row.push_back(new QStandardItem(formattedStatistics)); // Statistics

        model->appendRow(row);
    }

    stopCodeTimer("createModelFromData");

    return QVariant::fromValue(model);

}

QStringList SettingsAction::getSystemModeColor() {
    // Get the application palette
    QPalette palette = QApplication::palette();

    // Check the color of the window text to determine if the system is in dark mode or light mode
    // Assuming dark mode has lighter text (e.g., white) and light mode has darker text (e.g., black)
    if (palette.color(QPalette::WindowText).lightness() < 128) {
        // Light mode
        return { "#FFFFFF","#000000" }; // White
    }
    else {
        // Dark mode
        return { "#000000","#FFFFFF" }; // Black
    }
}

void SettingsAction::populatePointDataConcurrently(QString datasetId, const std::vector<float>& pointVector, int numPoints, int numDimensions, std::vector<QString> dimensionNames)
{
    QtConcurrent::run([this, datasetId, pointVector, numPoints, numDimensions, dimensionNames]() {
        auto pointDataset = mv::data().getDataset<Points>(datasetId);

        if (pointDataset.isValid())
        {
            pointDataset->setSelectionIndices({});
            if (!pointVector.empty() && numPoints > 0 && numDimensions > 0) {
                pointDataset->setData(pointVector.data(), numPoints, numDimensions);
                pointDataset->setDimensionNames(dimensionNames);
                mv::events().notifyDatasetDataChanged(pointDataset);
            }
        }
        });
}
void SettingsAction::enableActions()
{
    //_startComputationTriggerAction.setDisabled(false);
    _topNGenesFilter.setDisabled(false);
    _typeofTopNGenes.setDisabled(false);
    _clusterCountSortingType.setDisabled(false);
    _scatterplotReembedColorOption.setDisabled(false);
    _applyLogTransformation.setDisabled(false);
    _usePreComputedTSNE.setDisabled(false);
    _tsnePerplexity.setDisabled(false);
    _referenceTreeDataset.setDisabled(false);
    _mainPointsDataset.setDisabled(false);
    _embeddingDataset.setDisabled(false);
    _speciesNamesDataset.setDisabled(false);
    _clusterNamesDataset.setDisabled(false);
    _scatterplotEmbeddingPointsUMAPOption.setDisabled(false);
    _speciesExplorerInMap.setDisabled(false);
    _selectedSpeciesVals.setDisabled(false);
    _statusColorAction.setDisabled(false);
    _searchBox->setDisabled(false);
    enableDisableButtonsAutomatically();
    if (_statusColorAction.getString() == "C")
    {
        _startComputationTriggerAction.setDisabled(true);
    }
    else
    {
        _startComputationTriggerAction.setDisabled(false);
    }
}
void SettingsAction::disableActions()
{
    _statusColorAction.setString("R");
    _startComputationTriggerAction.setDisabled(true);
    _topNGenesFilter.setDisabled(true);
    _typeofTopNGenes.setDisabled(true);
    _clusterCountSortingType.setDisabled(true);
    _scatterplotReembedColorOption.setDisabled(true);
    _removeRowSelection.setDisabled(true);
    _revertRowSelectionChangesToInitial.setDisabled(true);
    _speciesExplorerInMapTrigger.setDisabled(true);
    _usePreComputedTSNE.setDisabled(true);
    _applyLogTransformation.setDisabled(true);
    _tsnePerplexity.setDisabled(true);
    _referenceTreeDataset.setDisabled(true);
    _mainPointsDataset.setDisabled(true);
    _embeddingDataset.setDisabled(true);
    _speciesNamesDataset.setDisabled(true);
    _clusterNamesDataset.setDisabled(true);
    _scatterplotEmbeddingPointsUMAPOption.setDisabled(true);
    _speciesExplorerInMap.setDisabled(true);
    _selectedSpeciesVals.setDisabled(true);
    _statusColorAction.setDisabled(true);
    _searchBox->setDisabled(true);
}

void SettingsAction::enableDisableButtonsAutomatically()
{

    bool optionsActionHasOptions = !_speciesExplorerInMap.getOptions().isEmpty();
    bool stringActionHasOptions = !_selectedSpeciesVals.getString().isEmpty();

    bool bothListsEqual = false;
    if (optionsActionHasOptions && stringActionHasOptions) {
        QStringList temp = _selectedSpeciesVals.getString().split(" @%$,$%@ ");
        QStringList species = _speciesExplorerInMap.getSelectedOptions();

        std::sort(temp.begin(), temp.end());
        std::sort(species.begin(), species.end());
        bothListsEqual = (temp == species);
    }

    if (!stringActionHasOptions)
    {
        _speciesExplorerInMapTrigger.setDisabled(true);
        _revertRowSelectionChangesToInitial.setDisabled(true);
    }
    else
    {
        if (!optionsActionHasOptions)
        {
            _speciesExplorerInMapTrigger.setDisabled(true);
            _revertRowSelectionChangesToInitial.setDisabled(false);
        }
        else
        {
            if (bothListsEqual)
            {
                _speciesExplorerInMapTrigger.setDisabled(false);
                _revertRowSelectionChangesToInitial.setDisabled(true);
            }
            else
            {
                _speciesExplorerInMapTrigger.setDisabled(false);
                _revertRowSelectionChangesToInitial.setDisabled(false);
            }
        }
    }



}


void SettingsAction::populatePointData(QString& datasetId, std::vector<float>& pointVector, int& numPoints, int& numDimensions, std::vector<QString>& dimensionNames)
{
    auto pointDataset = mv::data().getDataset<Points>(datasetId);

    if (pointDataset.isValid())
    {
        pointDataset->setSelectionIndices({});
        if (pointVector.size() > 0 && numPoints > 0 && numDimensions > 0) {
            pointDataset->setData(pointVector.data(), numPoints, numDimensions);
            pointDataset->setDimensionNames(dimensionNames);
            mv::events().notifyDatasetDataChanged(pointDataset);
        }

    }
}
void SettingsAction::populateClusterData(QString& datasetId, std::map<QString, std::pair<QColor, std::vector<int>>>& clusterMap)
{

    auto colorDataset = mv::data().getDataset<Clusters>(datasetId);
    if (colorDataset.isValid())
    {
        for (const auto& pair : clusterMap)
        {
            QString clusterName = pair.first;
            std::pair<QColor, std::vector<int>> value = pair.second;
            QColor clusterColor = value.first;
            std::vector<std::uint32_t> clusterIndices(value.second.begin(), value.second.end());

            if (clusterIndices.size() > 0)
            {
                Cluster clusterValue;
                clusterValue.setName(clusterName);
                clusterValue.setColor(clusterColor);
                clusterValue.setIndices(clusterIndices);
                colorDataset->addCluster(clusterValue);
            }
        }

        mv::events().notifyDatasetDataChanged(colorDataset);
    }


}

void SettingsAction::removeSelectionTableRows(QStringList* selectedLeaves)
{
    //check if _selectionDetailsTable is valid
    if (_selectionDetailsTable == nullptr) {
        return;
    }



    QAbstractItemModel* model = _selectionDetailsTable->model();

    //check if model is valid
    if (model == nullptr) {
        return;
    }

    auto colorValues = getSystemModeColor();
    auto systemColor = colorValues[0];
    auto valuesColor = colorValues[1];

    // Iterate through all rows
    for (int row = 0; row < model->rowCount(); ++row) {
        QModelIndex index = model->index(row, 0); // Assuming species name is in column 0
        QString species = model->data(index).toString();

        // Check if the species is one of the selected species
        if (selectedLeaves->contains(species)) {
            for (int col = 1; col < model->columnCount(); ++col) {
                QModelIndex cellIndex = model->index(row, col);
                _selectionDetailsTable->model()->setData(cellIndex, QBrush(QColor("#00A2ED")), Qt::BackgroundRole);
                _selectionDetailsTable->model()->setData(cellIndex, QBrush(QColor("#000000")), Qt::ForegroundRole);
            }
        }
        else
        {
            //remove existing color from rows
            for (int col = 1; col < model->columnCount(); ++col) {
                QModelIndex cellIndex = model->index(row, col);
                _selectionDetailsTable->model()->setData(cellIndex, QBrush(QColor(systemColor)), Qt::BackgroundRole);
                _selectionDetailsTable->model()->setData(cellIndex, QBrush(QColor(valuesColor)), Qt::ForegroundRole);
            }
        }
    }

}

void SettingsAction::updateSelectedSpeciesCounts(QJsonObject& node, const std::map<QString, int>& speciesCountMap) {
    // Check if the "name" key exists in the current node
    if (node.contains("name")) {
        QString nodeName = node["name"].toString();
        auto it = speciesCountMap.find(nodeName);
        // If the "name" is found in the speciesExpressionMap, update "mean" if it exists or add "mean" if it doesn't exist
        if (it != speciesCountMap.end()) {
            node["cellCounts"] = it->second; // Use it->second to access the value in the map
        }
    }

    // If the node has "children", recursively update them as well
    if (node.contains("children")) {
        QJsonArray children = node["children"].toArray();
        for (int i = 0; i < children.size(); ++i) {
            QJsonObject child = children[i].toObject();
            updateSelectedSpeciesCounts(child, speciesCountMap); // Recursive call
            children[i] = child; // Update the modified object back into the array
        }
        node["children"] = children; // Update the modified array back into the parent JSON object
    }
}
/*
QString SettingsAction::createJsonTreeFromNewick(QString tree, std::vector<QString> leafnames, std::map <QString, Statistics> speciesMeanValues)
{
    int i = 0;
    std::string jsonString = "";
    std::stringstream jsonStream;
    std::string newick = tree.toStdString();
    while (i < newick.size()) {
        if (newick[i] == '(') {
            jsonStream << "{\n\"children\": [";
            i++;
        }
        else if (newick[i] == ',') {
            jsonStream << ",";
            i++;
        }
        else if (newick[i] == ')') {
            jsonStream << "],\n\"id\": 1,\n\"score\": 1,\n\"branchLength\": 1.0,\n\"width\": 1\n}";
            i++;
        }
        else if (newick[i] == ';') {
            break;
        }
        else {
            if (isdigit(newick[i])) {
                int skip = 1;
                std::string num = "";
                for (int j = i; j < newick.size(); j++) {
                    if (isdigit(newick[j])) {
                        continue;
                    }
                    else {
                        num = newick.substr(i, j - i);

                        skip = j - i;
                        break;
                    }
                }
                std::string species = leafnames[(std::stoi(num) - 1)].toStdString();
                //std::string meanValue = std::to_string(speciesMeanValues[QString::fromStdString(species)]);
                auto it = speciesMeanValues.find(QString::fromStdString(species));
                std::string meanValue;
                if (it != speciesMeanValues.end()) {
                    // Key found, use the corresponding value
                    meanValue = std::to_string(it->second.meanSelected);
                }
                else {
                    // Key not found, assign -1
                    meanValue = "-1";
                }

                jsonStream << "{\n\"color\": \"#000000\",\n\"hastrait\": true,\n\"iscollapsed\": false,\n\"branchLength\": 1.0,\n\"cellCounts\": " << 0 << ",\n\"mean\": "<< meanValue <<", \n\"name\": \"" << species << "\"\n}";
                i += skip;
            }
        }
    }

    jsonString = jsonStream.str();

    nlohmann::json json = nlohmann::json::parse(jsonString);
    std::string jsonStr = json.dump(4);
    //qDebug()<< "CrossSpeciesComparisonClusterRankPlugin::createJsonTree: jsonStr: " << QString::fromStdString(jsonStr);
    QString formattedTree = QString::fromStdString(jsonStr);


    return  formattedTree;
}
*/
/*
std::string SettingsAction::mergeToNewick(int* merge, int numOfLeaves) {
    std::vector<std::string> labels(numOfLeaves);
    for (int i = 0; i < numOfLeaves; ++i) {
        labels[i] = std::to_string(i + 1);
    }

    std::stack<std::string> stack;

    for (int i = 0; i < 2 * (numOfLeaves - 1); i += 2) {
        int left = merge[i];
        int right = merge[i + 1];

        std::string leftStr;
        if (left < 0) {
            leftStr = labels[-left - 1];
        }
        else {
            leftStr = stack.top();
            stack.pop();
        }

        std::string rightStr;
        if (right < 0) {
            rightStr = labels[-right - 1];
        }
        else {
            rightStr = stack.top();
            stack.pop();
        }

        std::string merged = "(" + leftStr + "," + rightStr + ")";
        stack.push(merged);
    }

    return stack.top() + ";";
}
*/
double* SettingsAction::condensedDistanceMatrix(const std::vector<float>& items) {
    size_t n = items.size();
    double* distmat = new double[(n * (n - 1)) / 2];
    size_t k = 0;

#pragma omp parallel for schedule(dynamic) collapse(2) private(k)
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            k = ((n * (n - 1)) / 2) - ((n - i) * (n - i - 1)) / 2 + j - i - 1;
            distmat[k] = std::abs(items[i] - items[j]);
        }
    }

    return distmat;
}

SettingsAction::Widget::Widget(QWidget* parent, SettingsAction* SettingsAction) :
    WidgetActionWidget(parent, SettingsAction)
{ }

SettingsAction::OptionSelectionAction::Widget::Widget(QWidget* parent, OptionSelectionAction* optionSelectionAction) :
    WidgetActionWidget(parent, optionSelectionAction)
{ }

inline SettingsAction::OptionSelectionAction::OptionSelectionAction(SettingsAction& SettingsAction) :
    GroupAction(nullptr, "CrossSpeciesComparisonGeneDetectPluginOptionSelectionAction"),
    _settingsAction(SettingsAction)
{
    setText("Options");
    setIcon(Application::getIconFont("FontAwesome").getIcon("wrench"));
    //addAction(&_settingsAction.getTableModelAction());
    //addAction(&_settingsAction.getSelectedGeneAction());
    //addAction(&_settingsAction.getSelectedRowIndexAction());
    //addAction(&_settingsAction.getFilteringTreeDatasetAction());
    //addAction(&_settingsAction.getOptionSelectionAction());
    //addAction(&_settingsAction.getStartComputationTriggerAction());
    //addAction(&_settingsAction.getReferenceTreeDatasetAction());
    //addAction(&_settingsAction.getMainPointsDataset());
    //addAction(&_settingsAction.getHierarchyTopClusterDataset());
    //addAction(&_settingsAction.getHierarchyMiddleClusterDataset());
    //addAction(&_settingsAction.getHierarchyBottomClusterDataset());
    //addAction(&_settingsAction.getSpeciesNamesDataset());
    //addAction(&_settingsAction.getSelectedClusterNames());


}


void SettingsAction::fromVariantMap(const QVariantMap& variantMap)
{
    WidgetAction::fromVariantMap(variantMap);

    _geneNamesConnection.fromParentVariantMap(variantMap);
    _createRowMultiSelectTree.fromParentVariantMap(variantMap);
    _listModel.fromParentVariantMap(variantMap);
    _selectedGene.fromParentVariantMap(variantMap);
    _mainPointsDataset.fromParentVariantMap(variantMap);
    _embeddingDataset.fromParentVariantMap(variantMap);
    _speciesNamesDataset.fromParentVariantMap(variantMap);
    _clusterNamesDataset.fromParentVariantMap(variantMap);
    _filteredGeneNamesVariant.fromParentVariantMap(variantMap);
    _topNGenesFilter.fromParentVariantMap(variantMap);
    _filteringEditTreeDataset.fromParentVariantMap(variantMap);
    _referenceTreeDataset.fromParentVariantMap(variantMap);
    _selectedRowIndex.fromParentVariantMap(variantMap);
    _performGeneTableTsneAction.fromParentVariantMap(variantMap);
    _tsnePerplexity.fromParentVariantMap(variantMap);
    _hiddenShowncolumns.fromParentVariantMap(variantMap);
    _speciesExplorerInMap.fromParentVariantMap(variantMap);
    _scatterplotReembedColorOption.fromParentVariantMap(variantMap);
    _scatterplotEmbeddingPointsUMAPOption.fromParentVariantMap(variantMap);
    _selectedSpeciesVals.fromParentVariantMap(variantMap);
    _removeRowSelection.fromParentVariantMap(variantMap);
    _revertRowSelectionChangesToInitial.fromParentVariantMap(variantMap);
    _speciesExplorerInMapTrigger.fromParentVariantMap(variantMap);
    _statusColorAction.fromParentVariantMap(variantMap);
    _typeofTopNGenes.fromParentVariantMap(variantMap);
    _clusterCountSortingType.fromParentVariantMap(variantMap);
    _usePreComputedTSNE.fromParentVariantMap(variantMap);
    _applyLogTransformation.fromParentVariantMap(variantMap);

}

QVariantMap SettingsAction::toVariantMap() const
{
    QVariantMap variantMap = WidgetAction::toVariantMap();

    _geneNamesConnection.insertIntoVariantMap(variantMap);
    _createRowMultiSelectTree.insertIntoVariantMap(variantMap);
    _listModel.insertIntoVariantMap(variantMap);
    _selectedGene.insertIntoVariantMap(variantMap);
    _mainPointsDataset.insertIntoVariantMap(variantMap);
    _embeddingDataset.insertIntoVariantMap(variantMap);
    _speciesNamesDataset.insertIntoVariantMap(variantMap);
    _clusterNamesDataset.insertIntoVariantMap(variantMap);
    _filteredGeneNamesVariant.insertIntoVariantMap(variantMap);
    _topNGenesFilter.insertIntoVariantMap(variantMap);
    _filteringEditTreeDataset.insertIntoVariantMap(variantMap);
    _referenceTreeDataset.insertIntoVariantMap(variantMap);
    _selectedRowIndex.insertIntoVariantMap(variantMap);
    _performGeneTableTsneAction.insertIntoVariantMap(variantMap);
    _tsnePerplexity.insertIntoVariantMap(variantMap);
    _hiddenShowncolumns.insertIntoVariantMap(variantMap);
    _speciesExplorerInMap.insertIntoVariantMap(variantMap);
    _scatterplotReembedColorOption.insertIntoVariantMap(variantMap);
    _scatterplotEmbeddingPointsUMAPOption.insertIntoVariantMap(variantMap);
    _selectedSpeciesVals.insertIntoVariantMap(variantMap);
    _removeRowSelection.insertIntoVariantMap(variantMap);
    _revertRowSelectionChangesToInitial.insertIntoVariantMap(variantMap);
    _speciesExplorerInMapTrigger.insertIntoVariantMap(variantMap);
    _statusColorAction.insertIntoVariantMap(variantMap);
    _typeofTopNGenes.insertIntoVariantMap(variantMap);
    _clusterCountSortingType.insertIntoVariantMap(variantMap);
    _usePreComputedTSNE.insertIntoVariantMap(variantMap);
    _applyLogTransformation.insertIntoVariantMap(variantMap);
    return variantMap;
}
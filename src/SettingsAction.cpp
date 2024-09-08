#include "SettingsAction.h"
#include "CrossSpeciesComparisonGeneDetectPlugin.h"
#include<string>  
#include <QFileDialog>
#include <QPageLayout>
#include <QWebEngineView>
#include <CrossSpeciesComparisonTreeData.h>
#include <numeric>   // for std::reduce
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
#include <thread>
#include <unordered_map>
#include <QApplication>
#include <QPalette>
#ifdef _WIN32
#include <execution>
#endif
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
        sortOrderMap[customOrder[i]] = static_cast<int>(i);
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




Stats combineStatisticsSingle(const StatisticsSingle& selected, const StatisticsSingle& nonSelected,const int countTopAbundance/*, const StatisticsSingle& allSelected*/) {
    Stats combinedStats;
    combinedStats.meanSelected = selected.meanVal;
    combinedStats.countSelected = selected.countVal;

    //combinedStats.meanAll = allSelected.meanVal;
    //combinedStats.countAll = allSelected.countVal;

    combinedStats.meanNonSelected = nonSelected.meanVal; 
    combinedStats.countNonSelected = nonSelected.countVal; 

    combinedStats.abundanceCountTop = countTopAbundance;

    return combinedStats;
}

StatisticsSingle calculateStatistics(const std::vector<float>& data) {
    const int count = static_cast<int>(data.size());
    if (count == 0) {
        return { 0.0f, 0 }; // Return early if data is empty to avoid division by zero
    }

    float sum = 0.0;
#ifdef _WIN32
    sum = std::reduce(std::execution::par, data.begin(), data.end(), 0.0f);
#else
    sum = std::reduce(data.begin(), data.end(), 0.0f);
#endif
    //sum= std::accumulate(v.begin(), v.end(), 0.0);
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
        //float sum= std::accumulate(v.begin(), v.end(), 0.0);

        
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
    return (it != vec.end()) ? static_cast<int>(std::distance(vec.begin(), it)) : -1;
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
    _bottomClusterNamesDataset(this, "Bottom Cluster Names"),
    _middleClusterNamesDataset(this, "Middle Cluster Names"),
    _topClusterNamesDataset(this, "Top Cluster Names"),
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
    _topHierarchyClusterNamesFrequencyInclusionList(this, "Top Hierarchy Cluster Names Frequency Inclusion List"),
    _middleHierarchyClusterNamesFrequencyInclusionList(this, "Middle Hierarchy Cluster Names Frequency Inclusion List"),
    _bottomHierarchyClusterNamesFrequencyInclusionList(this, "Bottom Hierarchy Cluster Names Frequency Inclusion List"),
    _speciesExplorerInMapTrigger(this, "Explore"),
    _applyLogTransformation(this, "Gene mapping log"),
    _clusterCountSortingType(this, "Cluster Count Sorting Type"),
    _currentCellSelectionClusterInfoLabel(nullptr),
    _performGeneTableTsnePerplexity(this, "Perform Gene Table TSNE Perplexity"),
    _performGeneTableTsneKnn(this, "Perform Gene Table TSNE Knn"),
    _performGeneTableTsneDistance(this, "Perform Gene Table TSNE Distance"),
    _performGeneTableTsneTrigger(this, "Perform Gene Table TSNE Trigger"),
    _clusterOrderHierarchy(this, "Cluster Order Hierarchy"),
    _toggleScatterplotSelection(this, "Show Scatterplot Selection"),
    _mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker(this, "Map For Hierarchy Items Change Method Stop For Project Load Blocker")
{
    _mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker.setChecked(true);
    setSerializationName("CSCGDV:CrossSpeciesComparison Gene Detect Plugin Settings");
    _statusBarActionWidget = new QStatusBar();


    _popupMessage = new QMessageBox();
    _popupMessage->setIcon(QMessageBox::Information);
    _popupMessage->setWindowTitle("Computation in Progress");
    _popupMessage->setText("Data Precomputation in Progress");
    _popupMessage->setInformativeText("The system is currently precomputing essential data to enhance your interactive exploration experience. This process may take some time based on your data size and processor. The popup will close automatically once initialization is complete. Thank you for using Cytosplore EvoViewer.");
    _popupMessage->setStandardButtons(QMessageBox::NoButton);
    _popupMessage->setModal(true);


    /*
    _popupMessage = new QMessageBox();
    _popupMessage->setIcon(QMessageBox::Information);
    _popupMessage->setWindowTitle("Computation in Progress");
    _popupMessage->setText(
        "<div style='text-align: center;'>"
        "<h2 style='color: #333;'>Data Precomputation in Progress</h2>"
        "</div>"
    );
    _popupMessage->setInformativeText(
        "<div style='text-align: center;'>"
        "The system is currently precomputing essential data to enhance your interactive exploration experience. "
        "This process may take some time based on your data size and processor. The popup will close automatically "
        "once initialization is complete. Thank you for using Cytosplore EvoViewer."
        "</div>"
    );
    _popupMessage->setStandardButtons(QMessageBox::NoButton);
    _popupMessage->setModal(true);
    */

    _searchBox = new CustomLineEdit();
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
    _searchBox->setClearButtonEnabled(false);
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
    _geneTableView->setAlternatingRowColors(false);
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

    // removeDatasets(groupIDDeletion);
     //removeDatasets(groupId1);
     //removeDatasets(groupId2);
     //removeDatasets(groupId3);

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
    _mainPointsDataset.setSerializationName("CSCGDV:Main Points Dataset");
    _embeddingDataset.setSerializationName("CSCGDV:Embedding Dataset");
    _speciesNamesDataset.setSerializationName("CSCGDV:Species Names Dataset");
    _bottomClusterNamesDataset.setSerializationName("CSCGDV:Cluster Names Dataset");
    _middleClusterNamesDataset.setSerializationName("CSCGDV:Middle Cluster Names Dataset");
    _topClusterNamesDataset.setSerializationName("CSCGDV:Top Cluster Names Dataset");
    _filteredGeneNamesVariant.setSerializationName("CSCGDV:Filtered Gene Names");
    _topNGenesFilter.setSerializationName("CSCGDV:Top N Genes Filter");
    _filteringEditTreeDataset.setSerializationName("CSCGDV:Filtering Tree Dataset");
    _referenceTreeDataset.setSerializationName("CSCGDV:Reference Tree Dataset");
    _selectedRowIndex.setSerializationName("CSCGDV:Selected Row Index");
    _geneNamesConnection.setSerializationName("CSCGDV:Gene Names Connection");
    _selectedSpeciesVals.setSerializationName("CSCGDV:Selected Species Vals");
    _clusterOrderHierarchy.setSerializationName("CSCGDV:Cluster Order Hierarchy");
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
    _performGeneTableTsnePerplexity.setSerializationName("CSCGDV:Gene Table TSNE Perplexity");
    _performGeneTableTsnePerplexity.setMinimum(1);
    _performGeneTableTsnePerplexity.setMaximum(50);
    _performGeneTableTsnePerplexity.setValue(15);
    _performGeneTableTsneKnn.setSerializationName("CSCGDV:Gene Table TSNE Knn");
    _performGeneTableTsneKnn.initialize({ "FLANN","HNSW","ANNOY" }, "ANNOY");
    _performGeneTableTsneDistance.setSerializationName("CSCGDV:Gene Table TSNE Distance");
    _performGeneTableTsneDistance.initialize({ "Euclidean","Cosine","Inner Product","Manhattan","Hamming","Dot" }, "Dot");
    _performGeneTableTsneTrigger.setSerializationName("CSCGDV:Gene Table TSNE Trigger");
    _performGeneTableTsneTrigger.setDisabled(true);
    _clusterOrderHierarchy.setString("");
    _tsnePerplexity.setSerializationName("CSCGDV:TSNE Perplexity");
    _tsnePerplexity.setMinimum(1);
    _tsnePerplexity.setMaximum(50);
    _tsnePerplexity.setValue(30);
    _usePreComputedTSNE.setSerializationName("CSCGDV:Use Precomputed TSNE");
    _usePreComputedTSNE.setChecked(true);
    _applyLogTransformation.setChecked(false);
    _toggleScatterplotSelection.setChecked(false);
    _performGeneTableTsneAction.setChecked(false);
    _hiddenShowncolumns.setSerializationName("CSCGDV:Hidden Shown Columns");
    _speciesExplorerInMap.setSerializationName("CSCGDV:Species Explorer In Map");
    _topHierarchyClusterNamesFrequencyInclusionList.setSerializationName("CSCGDV:Top Hierarchy Cluster Names Frequency Inclusion List");
    _middleHierarchyClusterNamesFrequencyInclusionList.setSerializationName("CSCGDV:Middle Hierarchy Cluster Names Frequency Inclusion List");
    _bottomHierarchyClusterNamesFrequencyInclusionList.setSerializationName("CSCGDV:Bottom Hierarchy Cluster Names Frequency Inclusion List");
    _scatterplotReembedColorOption.setSerializationName("CSCGDV:Scatterplot Reembedding Color Option");
    _scatterplotEmbeddingPointsUMAPOption.setSerializationName("CSCGDV:Scatterplot Embedding UMAP Points Option");
    _typeofTopNGenes.setSerializationName("CSCGDV:Type of Top N Genes");
    _clusterCountSortingType.setSerializationName("CSCGDV:Cluster Count Sorting Type");
    _applyLogTransformation.setSerializationName("CSCGDV:Apply Log Transformation");
    _createRowMultiSelectTree.setDisabled(true);
    _selectedRowIndex.setDisabled(true);
    _selectedRowIndex.setString("");
    _scatterplotReembedColorOption.initialize({ "Species","Cluster","Expression" }, "Species");
    _typeofTopNGenes.initialize({ "Absolute","Negative","Positive" }, "Positive");
    _clusterCountSortingType.initialize({ "Count","Name","Hierarchy View" }, "Count");
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
    _bottomClusterNamesDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == ClusterType;
        });
    _middleClusterNamesDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == ClusterType;
        });
    _topClusterNamesDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
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

    const int delayMs = 500; // Delay in milliseconds

    QTimer* bottomTimer = new QTimer(this);
    bottomTimer->setSingleShot(true);
    const auto updateBottomHierarchyClusterNamesFrequencyInclusionList = [this, bottomTimer]() -> void
        {
            bottomTimer->start(delayMs);
        };
    connect(bottomTimer, &QTimer::timeout, this, [this]() {

        {
            computeFrequencyMapForHierarchyItemsChange("bottom");
            _statusColorAction.setString("M");
        }
        });
    connect(&_bottomHierarchyClusterNamesFrequencyInclusionList, &OptionsAction::selectedOptionsChanged, this, updateBottomHierarchyClusterNamesFrequencyInclusionList);

    QTimer* middleTimer = new QTimer(this);
    middleTimer->setSingleShot(true);
    const auto updateMiddleHierarchyClusterNamesFrequencyInclusionList = [this, middleTimer]() -> void
        {
            middleTimer->start(delayMs);
        };
    connect(middleTimer, &QTimer::timeout, this, [this]() {

        {
            computeFrequencyMapForHierarchyItemsChange("middle");
            _statusColorAction.setString("M");
        }

        });
    connect(&_middleHierarchyClusterNamesFrequencyInclusionList, &OptionsAction::selectedOptionsChanged, this, updateMiddleHierarchyClusterNamesFrequencyInclusionList);

    QTimer* topTimer = new QTimer(this);
    topTimer->setSingleShot(true);
    const auto updateTopHierarchyClusterNamesFrequencyInclusionList = [this, topTimer]() -> void
        {
            topTimer->start(delayMs);
        };
    connect(topTimer, &QTimer::timeout, this, [this]() {

        {
            computeFrequencyMapForHierarchyItemsChange("top");
            _statusColorAction.setString("M");
        }

        });
    connect(&_topHierarchyClusterNamesFrequencyInclusionList, &OptionsAction::selectedOptionsChanged, this, updateTopHierarchyClusterNamesFrequencyInclusionList);

    const auto updateGeneFilteringTrigger = [this]() -> void
        {
            disableActions();
            _pauseStatusUpdates = true;
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
            _totalGeneList.clear();
            auto fullDataset = mv::data().getDataset<Points>(_mainPointsDataset.getCurrentDataset().getDatasetId());
            auto dimensions = fullDataset->getNumDimensions();
            _totalGeneList = fullDataset->getDimensionNames();
            if (dimensions > 0) {
                _topNGenesFilter.setMinimum(1);
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
            const auto mainSelectionChanged = [this]() -> void {
                _toggleScatterplotSelection.setChecked(true);
                };
            connect(&fullDataset, &Dataset<Points>::dataSelectionChanged, this, mainSelectionChanged);

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
    const auto updateTopHierarchyDatasetChanged = [this]() -> void {
       
        if (_topClusterNamesDataset.getCurrentDataset().isValid())
        {
            auto clusterFullDataset = mv::data().getDataset<Clusters>(_topClusterNamesDataset.getCurrentDataset().getDatasetId());
            auto clusters = clusterFullDataset->getClusters();
            QStringList clusterNames = {};
            if (!clusters.isEmpty())
            {
                for (auto cluster : clusters)
                {
                    clusterNames.append(cluster.getName());
                }
            }
            _topHierarchyClusterNamesFrequencyInclusionList.setOptions(clusterNames);
            QString removalString = "Non-Neuronal";
            //if removal string present remove it
            if (clusterNames.contains(removalString))
            {
                clusterNames.removeAll(removalString);
               
            }
            _topHierarchyClusterNamesFrequencyInclusionList.setSelectedOptions(clusterNames);
        }
        else
        {
            _topHierarchyClusterNamesFrequencyInclusionList.setOptions({});
        }
        };

    connect(&_topClusterNamesDataset, &DatasetPickerAction::currentIndexChanged, this, updateTopHierarchyDatasetChanged);

    const auto updateMiddleHierarchyDatasetChanged = [this]() -> void {
        if (_middleClusterNamesDataset.getCurrentDataset().isValid())
        {
            auto clusterFullDataset = mv::data().getDataset<Clusters>(_middleClusterNamesDataset.getCurrentDataset().getDatasetId());
            auto clusters = clusterFullDataset->getClusters();
            QStringList clusterNames = {};
            if (!clusters.isEmpty())
            {
                for (auto cluster : clusters)
                {
                    clusterNames.append(cluster.getName());
                }
            }
            _middleHierarchyClusterNamesFrequencyInclusionList.setOptions(clusterNames);
            //QStringList removalStringList = { "Oligo","Astro", "OPC", "Micro-PVM, "Endo", "VLMC" };
            QStringList removalStringList = { "Oligo", "Astro", "OPC", "Micro-PVM", "Endo", "VLMC" };

            //if removal string present remove it
            for (auto removalString : removalStringList)
            {
                if (clusterNames.contains(removalString))
            {
                clusterNames.removeAll(removalString);
                
            }
                }
            
            _middleHierarchyClusterNamesFrequencyInclusionList.setSelectedOptions(clusterNames);
        }
        else
        {
            _middleHierarchyClusterNamesFrequencyInclusionList.setOptions({});
        }
        };

    connect(&_middleClusterNamesDataset, &DatasetPickerAction::currentIndexChanged, this, updateMiddleHierarchyDatasetChanged);

    const auto updateBottomHierarchyDatasetChanged = [this]() -> void {
        if (_bottomClusterNamesDataset.getCurrentDataset().isValid())
        {
            auto clusterFullDataset = mv::data().getDataset<Clusters>(_bottomClusterNamesDataset.getCurrentDataset().getDatasetId());
            auto clusters = clusterFullDataset->getClusters();
            QStringList clusterNames = {};
            if (!clusters.isEmpty())
            {
                for (auto cluster : clusters)
                {
                    clusterNames.append(cluster.getName());
                }
            }
            _bottomHierarchyClusterNamesFrequencyInclusionList.setOptions(clusterNames);
            QStringList removalStringList = { "Oligo_2","Oligo_1","Astro_2","Astro_1", "OPC", "Microglia/PVM", "Endo", "VLMC","LMC"};

            //if removal string present remove it
            for (auto removalString : removalStringList)
            {
                if (clusterNames.contains(removalString))
                {
                    clusterNames.removeAll(removalString);

                }
            }
                _bottomHierarchyClusterNamesFrequencyInclusionList.setSelectedOptions(clusterNames);
            
        }
        else
        {
            _bottomHierarchyClusterNamesFrequencyInclusionList.setOptions({});
        }
        };

    connect(&_bottomClusterNamesDataset, &DatasetPickerAction::currentIndexChanged, this, updateBottomHierarchyDatasetChanged);

    const auto updateScatterplotColor = [this]() -> void {
        auto selectedColorType= _scatterplotReembedColorOption.getCurrentText();
        if (selectedColorType != "")
        {
            auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
            mv::gui::DatasetPickerAction* colorDatasetPickerAction;
            mv::gui::DatasetPickerAction* pointDatasetPickerAction;
            mv::gui::ViewPluginSamplerAction* samplerActionAction;
            if (scatterplotViewFactory) {
                for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                    if (plugin->getGuiName() == "Scatterplot Cell Selection Overview") {
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
                                        if (_bottomClusterNamesDataset.getCurrentDataset().isValid())
                                        {
                                            colorDatasetPickerAction->setCurrentText("");
                                            colorDatasetPickerAction->setCurrentDataset(_bottomClusterNamesDataset.getCurrentDataset());
                                        }
                                    }
                                    else if (selectedColorType == "Species")
                                    {
                                        if (_speciesNamesDataset.getCurrentDataset().isValid())
                                        {
                                            colorDatasetPickerAction->setCurrentText("");
                                            colorDatasetPickerAction->setCurrentDataset(_speciesNamesDataset.getCurrentDataset());
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
                            
                            samplerActionAction = plugin->findChildByPath<mv::gui::ViewPluginSamplerAction>("Sampler");

                            if (samplerActionAction)
                            {
                                samplerActionAction->setTooltipGeneratorFunction([this](const ViewPluginSamplerAction::SampleContext& toolTipContext) -> QString {
                                    QString clusterDatasetId = _speciesNamesDataset.getCurrentDataset().getDatasetId();
                                    return generateTooltip(toolTipContext, clusterDatasetId,true, "GlobalPointIndices");
                                    });
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
        if (_pauseStatusUpdates)
        {
            return;
        }
        auto string = _statusColorAction.getString();
        QString labelText = "";
        QString backgroundColor = "none";
        if (string == "C")
        {
            _startComputationTriggerAction.setDisabled(true);
            if (_popupMessage->isVisible())
            {
                _popupMessage->hide();
            }
        }
        else
        {
            _startComputationTriggerAction.setDisabled(false);
        }
        if (string == "M")
        {
            _removeRowSelection.trigger();
        }

        if (string == "C") {
            labelText = "Updated";
            backgroundColor = "#28a745"; // Green
            
        }
        else if (string == "M") {
            labelText = "Outdated";
            backgroundColor = "#ffc107"; // Gold

        }
        else if (string == "E") {
            labelText ="Error";
            backgroundColor = "#dc3545"; // Red
        }
        else if (string == "R")
        {
            labelText = "Processing";
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


    const auto updateClusterOrderHierarchy = [this]() -> void {
        
        _customOrderClustersFromHierarchy.clear();
        if (_clusterOrderHierarchy.getString() != "")
        {
            QStringList clusterOrderHierarchyList = _clusterOrderHierarchy.getString().split(" @%$,$%@ ");
            for (auto clusterOrderHierarchyItem : clusterOrderHierarchyList)
            {

                _customOrderClustersFromHierarchy.push_back(clusterOrderHierarchyItem);
            }
        }
        
        


        };
    connect(&_clusterOrderHierarchy, &StringAction::changed, this, updateClusterOrderHierarchy);



    const auto updateApplyLogTransformation = [this]() -> void {
        _statusColorAction.setString("M");

        };
    connect(&_applyLogTransformation, &ToggleAction::toggled, this, updateApplyLogTransformation);
    const auto updateMapForHierarchyItemsChangeMethodStopForProjectLoadBlocker = [this]() -> void {
        if (!_mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker.isChecked())
        {
            _startComputationTriggerAction.setDisabled(false);
           
            _popupMessage->show();
            QApplication::processEvents();
            QFuture<void> future1 = QtConcurrent::run([this]() { computeGeneMeanExpressionMap(); });
            //QFuture<void> future3 = QtConcurrent::run([this]() { triggerTrippleHierarchyFrequencyChange(); });
            //QFuture<void> future4 = QtConcurrent::run([this]() { computeFrequencyMapForHierarchyItemsChange("top"); });
            //QFuture<void> future5 = QtConcurrent::run([this]() { computeFrequencyMapForHierarchyItemsChange("middle"); });
            //QFuture<void> future6 = QtConcurrent::run([this]() { computeFrequencyMapForHierarchyItemsChange("bottom"); });
            computeFrequencyMapForHierarchyItemsChange("top");
            future1.waitForFinished();
            precomputeTreesFromHierarchy();
            //computeFrequencyMapForHierarchyItemsChange("middle");
            //computeFrequencyMapForHierarchyItemsChange("bottom");
            

            _startComputationTriggerAction.trigger();
        }
        else
        {
            _startComputationTriggerAction.setDisabled(true);
        }

        };
    connect(&_mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker, &ToggleAction::toggled, this, updateMapForHierarchyItemsChangeMethodStopForProjectLoadBlocker);
        const auto updateToggleScatterplotSelection = [this]() -> void {
        
            auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
            mv::gui::DecimalAction* overlayopacityAction;

            if (scatterplotViewFactory) {
                for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                    if (plugin->getGuiName() == "Scatterplot Embedding View") {
                        overlayopacityAction = dynamic_cast<DecimalAction*>(plugin->findChildByPath("Settings/Selection/Opacity"));
                        if (overlayopacityAction)
                        {
                            //qDebug() << "Overlay opacity action found";
                            if (_toggleScatterplotSelection.isChecked())
                            {
                                overlayopacityAction->setValue(100.0);
                            }
                            else
                            {
                                overlayopacityAction->setValue(0.0);
                            }   
                        }
                    }
                }
            }

        };
    connect(&_toggleScatterplotSelection, &ToggleAction::toggled, this, updateToggleScatterplotSelection);


    const auto recomputeGeneTableTSNE = [this]() -> void {
        if (_selectedPointsTSNEDatasetForGeneTable.isValid())
        {
            
            auto runningAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDatasetForGeneTable->findChildByPath("TSNE/TsneComputationAction/Running"));

            if (runningAction)
            {

                if (runningAction->isChecked())
                {
                    auto stopAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDatasetForGeneTable->findChildByPath("TSNE/TsneComputationAction/Stop"));
                    if (stopAction)
                    {
                        stopAction->trigger();
                        std::this_thread::sleep_for(std::chrono::seconds(5));
                    }
                }

            }
            

            auto startAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDatasetForGeneTable->findChildByPath("TSNE/TsneComputationAction/Start"));
            if (startAction) {

                startAction->trigger();
            }

        }

        };
    connect(&_performGeneTableTsneTrigger, &TriggerAction::triggered, this, recomputeGeneTableTSNE);

    const auto updateGeneTableTSNECheck = [this]() -> void {
        _statusColorAction.setString("M");

        };
    connect(&_performGeneTableTsneAction, &ToggleAction::toggled, this, updateGeneTableTSNECheck);
    const auto updateClusterCountSortingType = [this]() -> void {
        updateClusterInfoStatusBar();

        };
    connect(&_clusterCountSortingType, &OptionAction::currentIndexChanged, this, updateClusterCountSortingType);
    QTimer* debounceTimer = new QTimer(this);
    debounceTimer->setSingleShot(true);
    debounceTimer->setInterval(500); // 500 milliseconds wait time

    const auto debouncelambda = [this]() -> void { // Capture debounceTimer by 
        disableActions();
        findTopNGenesPerCluster();
        enableActions();
        };

    connect(debounceTimer, &QTimer::timeout, this, debouncelambda);

    const auto updateTopGenesSlider = [this, debounceTimer]() -> void { // Capture debounceTimer by reference
        //wait to see if any more updates are coming then call findTopNGenesPerCluster();
        // Restart the timer every time the value changes
        debounceTimer->start();
        };
    connect(&_topNGenesFilter, &IntegralAction::valueChanged, this, updateTopGenesSlider);


    _statusColorAction.setString("M");
}

void SettingsAction::triggerTrippleHierarchyFrequencyChange()
{
    if (_mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker.isChecked())
    {
        return;
    }
    _clusterSpeciesFrequencyMap.clear();
    auto startTimer = std::chrono::high_resolution_clock::now();
    qDebug() << "computeFrequencyMapForHierarchyItemsChange for all 3 levels Start";

    if (!_speciesNamesDataset.getCurrentDataset().isValid() || !_mainPointsDataset.getCurrentDataset().isValid() || !_topClusterNamesDataset.getCurrentDataset().isValid() || !_middleClusterNamesDataset.getCurrentDataset().isValid() || !_bottomClusterNamesDataset.getCurrentDataset().isValid()) {
        qDebug() << "Datasets are not valid";
        return;
    }

    auto speciesClusterDatasetFull = mv::data().getDataset<Clusters>(_speciesNamesDataset.getCurrentDataset().getDatasetId());
    auto mainPointDatasetFull = mv::data().getDataset<Points>(_mainPointsDataset.getCurrentDataset().getDatasetId());
    auto numOfPoints = mainPointDatasetFull->getNumPoints();
    std::vector<bool> topClusterNames(numOfPoints, true);
    std::vector<bool> middleClusterNames(numOfPoints, true);
    std::vector<bool> bottomClusterNames(numOfPoints, true);
    QStringList topInclusionList;
    QStringList middleInclusionList;
    QStringList bottomInclusionList;
    auto topClusterDataset = mv::data().getDataset<Clusters>(_topClusterNamesDataset.getCurrentDataset().getDatasetId());
    auto middleClusterDataset = mv::data().getDataset<Clusters>(_middleClusterNamesDataset.getCurrentDataset().getDatasetId());
    auto bottomClusterDataset = mv::data().getDataset<Clusters>(_bottomClusterNamesDataset.getCurrentDataset().getDatasetId());

    auto processTopClusters = [&]() {
        if (topClusterDataset.isValid())
        {
            for (const auto& cluster : topClusterDataset->getClusters())
            {
                if (!topInclusionList.contains(cluster.getName()))
                {
                    for (const auto& index : cluster.getIndices())
                    {
                        topClusterNames[index] = false;
                    }
                }
            }
        }
        };

    auto processMiddleClusters = [&]() {
        if (middleClusterDataset.isValid())
        {
            for (const auto& cluster : middleClusterDataset->getClusters())
            {
                if (!middleInclusionList.contains(cluster.getName()))
                {
                    for (const auto& index : cluster.getIndices())
                    {
                        middleClusterNames[index] = false;
                    }
                }
            }
        }
        };

    auto processBottomClusters = [&]() {
        if (bottomClusterDataset.isValid())
        {
            for (const auto& cluster : bottomClusterDataset->getClusters())
            {
                if (!bottomInclusionList.contains(cluster.getName()))
                {
                    for (const auto& index : cluster.getIndices())
                    {
                        bottomClusterNames[index] = false;
                    }
                }
            }
        }
        };

    // Run the three tasks in parallel
    QFuture<void> topFuture = QtConcurrent::run(processTopClusters);
    QFuture<void> middleFuture = QtConcurrent::run(processMiddleClusters);
    QFuture<void> bottomFuture = QtConcurrent::run(processBottomClusters);

    // Wait for all tasks to complete
    topFuture.waitForFinished();
    middleFuture.waitForFinished();
    bottomFuture.waitForFinished();

    if (speciesClusterDatasetFull.isValid() && mainPointDatasetFull.isValid())
    {
        auto speciesclusters = speciesClusterDatasetFull->getClusters();
        for (const auto& species : speciesclusters) {
            auto speciesIndices = species.getIndices();
            auto speciesName = species.getName();
            int topCount = std::count_if(speciesIndices.begin(), speciesIndices.end(), [&topClusterNames](int index) {
                return topClusterNames[index];
                });
            int middleCount = std::count_if(speciesIndices.begin(), speciesIndices.end(), [&middleClusterNames](int index) {
                return middleClusterNames[index];
                });
            int bottomCount = std::count_if(speciesIndices.begin(), speciesIndices.end(), [&bottomClusterNames](int index) {
                return bottomClusterNames[index];
                });

            _clusterSpeciesFrequencyMap[speciesName]["topCells"] = topCount;
            _clusterSpeciesFrequencyMap[speciesName]["middleCells"] = middleCount;
            _clusterSpeciesFrequencyMap[speciesName]["bottomCells"] = bottomCount;
        }
    }

    auto endTimer = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTimer - startTimer).count();
    qDebug() << "Time taken for computeFrequencyMapForHierarchyItemsChange for all 3 levels: " + QString::number(duration / 1000.0) + " s";
}

void SettingsAction::updateButtonTriggered()
{
    if (_mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker.isChecked())
    {
        qDebug() << "Map for hierarchy items change method stop for project load blocker is checked";
        return;
    }

    try {
       // _startComputationTriggerAction.setDisabled(true);
        startCodeTimer("UpdateGeneFilteringTrigger");
        //startCodeTimer("Part1");

        int groupIDDeletion = 10;
        int groupID1 = 10 * 2;
        int groupID2 = 10 * 3;
        removeDatasets(groupIDDeletion);
        auto pointsDataset = _mainPointsDataset.getCurrentDataset();
        auto embeddingDataset = _embeddingDataset.getCurrentDataset();
        auto speciesDataset = _speciesNamesDataset.getCurrentDataset();
        auto clusterDataset = _bottomClusterNamesDataset.getCurrentDataset();
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

                    //if (_selectedPointsTSNEDataset.isValid())
                    //{
                        //auto datasetIDLowRem = _selectedPointsTSNEDataset.getDatasetId();
                        //mv::events().notifyDatasetAboutToBeRemoved(_selectedPointsTSNEDataset);
                        //mv::data().removeDataset(_selectedPointsTSNEDataset);
                        //mv::events().notifyDatasetRemoved(datasetIDLowRem, PointType);
                    //}
                    
                   // _selectedPointsDataset = Dataset<Points>();
                    //_selectedPointsEmbeddingDataset = Dataset<Points>();
                    //startCodeTimer("Part6.1");
                    /*if (!_selectedPointsDataset.isValid())
                    {
                        _selectedPointsDataset = mv::data().createDataset("Points", "SelectedPointsDataset");
                        _selectedPointsDataset->setGroupIndex(10);
                        mv::events().notifyDatasetAdded(_selectedPointsDataset);

                    }*/
                    

                    pointsDatasetRaw->setSelectionIndices(_selectedIndicesFromStorage);
                    _selectedPointsDataset = pointsDatasetRaw->createSubsetFromSelection("SelectedPointsDataset");
                    _selectedPointsDataset->setGroupIndex(groupIDDeletion);

                    if (!_tsneDatasetExpressionColors.isValid())
                    {
                        _tsneDatasetExpressionColors = mv::data().createDataset("Points", "TSNEDatasetExpressionColors", _selectedPointsDataset);
                        _tsneDatasetExpressionColors->setGroupIndex(groupIDDeletion);
                        mv::events().notifyDatasetAdded(_tsneDatasetExpressionColors);

                    }

                    embeddingDatasetRaw->setSelectionIndices(_selectedIndicesFromStorage);
                    _selectedPointsEmbeddingDataset = embeddingDatasetRaw->createSubsetFromSelection("TSNEDataset", _selectedPointsDataset);
                    _selectedPointsEmbeddingDataset->setGroupIndex(groupIDDeletion);


                    if (!_tsneDatasetSpeciesColors.isValid())
                    {
                        _tsneDatasetSpeciesColors = mv::data().createDataset("Cluster", "TSNEDatasetSpeciesColors", _selectedPointsDataset);
                        _tsneDatasetSpeciesColors->setGroupIndex(groupIDDeletion);
                        mv::events().notifyDatasetAdded(_tsneDatasetSpeciesColors);
                    }

                    if (!_tsneDatasetClusterColors.isValid())
                    {
                        _tsneDatasetClusterColors = mv::data().createDataset("Cluster", "TSNEDatasetClusterColors", _selectedPointsDataset);
                        _tsneDatasetClusterColors->setGroupIndex(groupIDDeletion);
                        mv::events().notifyDatasetAdded(_tsneDatasetClusterColors);
                    }


                    if (!_filteredUMAPDatasetPoints.isValid())
                    {
                        _filteredUMAPDatasetPoints = mv::data().createDataset("Points", "Filtered UMAP Dataset Points");
                        _filteredUMAPDatasetPoints->setGroupIndex(groupID1);
                        mv::events().notifyDatasetAdded(_filteredUMAPDatasetPoints);
                        if (!_filteredUMAPDatasetColors.isValid())
                        {
                            //need to delete

                        }
                        if (!_filteredUMAPDatasetClusters.isValid())
                        {
                            //need to delete

                        }
                        _filteredUMAPDatasetColors = mv::data().createDataset("Points", "Filtered UMAP Dataset Colors", _filteredUMAPDatasetPoints);
                        _filteredUMAPDatasetColors->setGroupIndex(groupID1);
                        mv::events().notifyDatasetAdded(_filteredUMAPDatasetColors);

                        _filteredUMAPDatasetClusters = mv::data().createDataset("Cluster", "Filtered UMAP Dataset Clusters", _filteredUMAPDatasetPoints);
                        _filteredUMAPDatasetClusters->setGroupIndex(groupID1);
                        mv::events().notifyDatasetAdded(_filteredUMAPDatasetClusters);

                    }

                    
                    
                    /*if (!_selectedPointsEmbeddingDataset.isValid())
                    {
                        _selectedPointsEmbeddingDataset = mv::data().createDataset("Points", "TSNEDataset", _selectedPointsDataset);
                        _selectedPointsEmbeddingDataset->setGroupIndex(10);
                        mv::events().notifyDatasetAdded(_selectedPointsEmbeddingDataset);

                    }*/


                    if (!_geneSimilarityPoints.isValid())
                    {
                        _geneSimilarityPoints = mv::data().createDataset("Points", "GeneSimilarityPoints");
                        _geneSimilarityPoints->setGroupIndex(groupID2);
                        mv::events().notifyDatasetAdded(_geneSimilarityPoints);
                    }
                    if (!_geneSimilarityClusterColoring.isValid())
                    {
                        _geneSimilarityClusterColoring = mv::data().createDataset("Cluster", "GeneSimilarityClusterColoring", _geneSimilarityPoints);
                        _geneSimilarityClusterColoring->setGroupIndex(groupID2);
                        mv::events().notifyDatasetAdded(_geneSimilarityClusterColoring);

                    }
                    //_geneSimilarityClusters.clear();
                    //stopCodeTimer("Part6.1");
                    if (_selectedPointsDataset.isValid() && _selectedPointsEmbeddingDataset.isValid() && _tsneDatasetSpeciesColors.isValid() && _tsneDatasetClusterColors.isValid() && _geneSimilarityPoints.isValid() && _geneSimilarityClusterColoring.isValid())
                    {
                        //startCodeTimer("Part6.2");
                        //_tsneDatasetSpeciesColors->getClusters() = QVector<Cluster>();
                        //events().notifyDatasetDataChanged(_tsneDatasetSpeciesColors);
                        //_tsneDatasetClusterColors->getClusters() = QVector<Cluster>();
                        //events().notifyDatasetDataChanged(_tsneDatasetClusterColors);
                        _geneSimilarityClusterColoring->getClusters() = QVector<Cluster>();
                        events().notifyDatasetDataChanged(_geneSimilarityClusterColoring);
                        //stopCodeTimer("Part6.2");
                        //startCodeTimer("Part7");
                        //startCodeTimer("Part7.1");
                        int selectedIndicesFromStorageSize = static_cast<int>(_selectedIndicesFromStorage.size());
                        int pointsDatasetColumnsSize = static_cast<int>(pointsDatasetallColumnIndices.size());
                        int embeddingDatasetColumnsSize = static_cast<int>(embeddingDatasetColumnIndices.size());
                        //QString datasetIdEmb = _selectedPointsDataset->getId();
                        //QString datasetId = _selectedPointsEmbeddingDataset->getId();
                        int dimofDatasetExp = 1;
                        std::vector<QString> dimensionNamesExp = { "Expression" };
                        QString datasetIdExp = _tsneDatasetExpressionColors->getId();
                        //stopCodeTimer("Part7.1");
                        //startCodeTimer("Part7.2");

                        // Define result containers outside the lambda functions to ensure they are accessible later
                        //std::vector<float> resultContainerForSelectedPoints(selectedIndicesFromStorageSize * pointsDatasetColumnsSize);
                        //std::vector<float> resultContainerForSelectedEmbeddingPoints(selectedIndicesFromStorageSize * embeddingDatasetColumnsSize);
                        std::vector<float> resultContainerColorPoints(selectedIndicesFromStorageSize, -1.0f);

                        //first thread start
                        //auto future1 = std::async(std::launch::async, [&]() {
                            //pointsDatasetRaw->populateDataForDimensions(resultContainerForSelectedPoints, pointsDatasetallColumnIndices, _selectedIndicesFromStorage);
                           // });

                        //second thread start
                       // auto future2 = std::async(std::launch::async, [&]() {
                            //embeddingDatasetRaw->populateDataForDimensions(resultContainerForSelectedEmbeddingPoints, embeddingDatasetColumnIndices, _selectedIndicesFromStorage);
                           // });


                        // Wait for all futures to complete before proceeding
                        //future1.wait();
                        //future2.wait();


                        //startCodeTimer("Part7.2.1");
                        //needs to wait for future1 finish only
                        //populatePointData(datasetIdEmb, resultContainerForSelectedPoints, selectedIndicesFromStorageSize, pointsDatasetColumnsSize, pointsDatasetallColumnNameList);
                        //stopCodeTimer("Part7.2.1");

                        //startCodeTimer("Part7.2.2");
                        //needs to wait for future2 finish only
                        //populatePointData(datasetId, resultContainerForSelectedEmbeddingPoints, selectedIndicesFromStorageSize, embeddingDatasetColumnsSize, embeddingDatasetallColumnNameList);
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
                        //stopCodeTimer("Part8");
                        //startCodeTimer("Part9");
                        if (!_performGeneTableTsneAction.isChecked())
                        {

                        
                        mv::plugin::AnalysisPlugin* analysisPlugin;
                        bool usePreTSNE = _usePreComputedTSNE.isChecked();

                        auto scatterplotModificationsLowDimUMAP = [this]() {
                            if (_selectedPointsTSNEDataset.isValid()) {
                                auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
                                mv::gui::DatasetPickerAction* colorDatasetPickerAction;
                                mv::gui::DatasetPickerAction* pointDatasetPickerAction;
                                mv::gui::ViewPluginSamplerAction* samplerActionAction;
                                if (scatterplotViewFactory) {
                                    for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                                        if (plugin->getGuiName() == "Scatterplot Cell Selection Overview") {
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
                                                            if (_bottomClusterNamesDataset.getCurrentDataset().isValid())
                                                            {
                                                                colorDatasetPickerAction->setCurrentDataset(_bottomClusterNamesDataset.getCurrentDataset());
                                                            }
                                                        }
                                                        else if (selectedColorType == "Species")
                                                        {
                                                            if (_speciesNamesDataset.getCurrentDataset().isValid())
                                                            {
                                                                colorDatasetPickerAction->setCurrentDataset(_speciesNamesDataset.getCurrentDataset());
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
                                                
                                                samplerActionAction = plugin->findChildByPath<mv::gui::ViewPluginSamplerAction>("Sampler");

                                                if (samplerActionAction)
                                                {
                                                    samplerActionAction->setTooltipGeneratorFunction([this](const ViewPluginSamplerAction::SampleContext& toolTipContext) -> QString {
                                                        QString clusterDatasetId = _speciesNamesDataset.getCurrentDataset().getDatasetId();
                                                        return generateTooltip(toolTipContext, clusterDatasetId,true, "GlobalPointIndices");
                                                        });
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
                            _selectedPointsTSNEDataset->setGroupIndex(groupIDDeletion);
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
                                _selectedPointsTSNEDataset->setGroupIndex(groupIDDeletion);
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
                                const auto& nonSelectionDetails = _clusterGeneMeanExpressionMap[speciesName][geneName]["allCells"];
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


                                int topHierarchyCountValue = 0;
                                if (_clusterSpeciesFrequencyMap.find(speciesName) != _clusterSpeciesFrequencyMap.end())
                                {
                                    topHierarchyCountValue = _clusterSpeciesFrequencyMap[speciesName]["topCells"];
                                }
                                float topHierarchyFrequencyValue = 0.0;
                                if (topHierarchyCountValue != 0.0f) {
                                    topHierarchyFrequencyValue = static_cast<float>(calculateStatisticsShort.countVal) / topHierarchyCountValue;
                                }



                                localClusterNameToGeneNameToExpressionValue[geneName] = combineStatisticsSingle(calculateStatisticsShort, calculateStatisticsNot, topHierarchyCountValue);
                            }

                            // Merge results in a thread-safe manner
                            QMutexLocker locker(&clusterNameToGeneNameToExpressionValueMutex);
                            for (const auto& pair : localClusterNameToGeneNameToExpressionValue) {
                                _clusterNameToGeneNameToExpressionValue[speciesName][pair.first] = pair.second;
                                _selectedSpeciesCellCountMap[speciesName].selectedCellsCount = pair.second.countSelected;
                                _selectedSpeciesCellCountMap[speciesName].nonSelectedCellsCount = pair.second.countNonSelected;
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
                    updateClusterInfoStatusBar();
                    /*
                    QLayoutItem* layoutItem;
                    while ((layoutItem = _selectedCellClusterInfoStatusBar->takeAt(0)) != nullptr) {
                        delete layoutItem->widget();
                        delete layoutItem;
                    }
                    if (_tsneDatasetClusterColors.isValid())
                    {

                        auto clusterValues = _tsneDatasetClusterColors->getClusters();
                        if (!clusterValues.empty())
                        {
                            //startCodeTimer("Part13");

                            //QLayoutItem* layoutItem;
                            //while ((layoutItem = _selectedCellClusterInfoStatusBar->takeAt(0)) != nullptr) {
                            //    delete layoutItem->widget();
                            //    delete layoutItem;
                            //}

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





                            //for (auto cluster : clusterValues) {
                            //    auto clusterName = cluster.getName();
                            //    auto clusterIndicesSize = cluster.getIndices().size();
                            //    auto clusterColor = cluster.getColor(); // Assuming getColor() returns a QColor

                            //    // Calculate luminance
                            //    qreal luminance = 0.299 * clusterColor.redF() + 0.587 * clusterColor.greenF() + 0.114 * clusterColor.blueF();

                            //    // Choose text color based on luminance
                            //    QString textColor = (luminance > 0.5) ? "black" : "white";

                            //    // Convert QColor to hex string for stylesheet
                            //    QString backgroundColor = clusterColor.name(QColor::HexArgb);

                            //    auto clusterLabel = new QLabel(QString("%1: %2").arg(clusterName).arg(clusterIndicesSize));
                            //    // Add text color and background color to clusterLabel with padding and border for better styling
                            //    clusterLabel->setStyleSheet(QString("QLabel { color: %1; background-color: %2; padding: 2px; border: 0.5px solid %3; }")
                            //        .arg(textColor).arg(backgroundColor).arg(textColor));
                            //    _selectedCellClusterInfoStatusBar->addWidget(clusterLabel);
                            //}


                        }

                    }
                    */
                    //the next line should only execute if all above are finished


                    //startCodeTimer("Part14");
                     findTopNGenesPerCluster();
                    //stopCodeTimer("Part14");

                    


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
        _statusColorAction.setString("E");
    }
    catch (...) {
        qDebug() << "An unknown exception occurred in coputation";
        _statusColorAction.setString("E");
    }
}

void SettingsAction::updateClusterInfoStatusBar()
{
    QLayoutItem* layoutItem;
    while ((layoutItem = _selectedCellClusterInfoStatusBar->takeAt(0)) != nullptr) {
        delete layoutItem->widget();
        delete layoutItem;
    }
    if (_tsneDatasetClusterColors.isValid() && _bottomClusterNamesDataset.getCurrentDataset().isValid())
    {
        auto clusterDatasetName = _bottomClusterNamesDataset.getCurrentDataset()->getGuiName();
        auto clusterValues = _tsneDatasetClusterColors->getClusters();
        if (!clusterValues.empty())
        {
            //startCodeTimer("Part13");

            /*QLayoutItem* layoutItem;
            while ((layoutItem = _selectedCellClusterInfoStatusBar->takeAt(0)) != nullptr) {
                delete layoutItem->widget();
                delete layoutItem;
            }*/

            // Create a description label
            auto descriptionLabel = new QLabel("Cell counts per " + clusterDatasetName + ", sorted by " + _clusterCountSortingType.getCurrentText() + ":");

            // Optionally, set a stylesheet for the description label for styling
            descriptionLabel->setStyleSheet("QLabel { font-weight: bold; padding: 2px; }");
            // Add the description label to the layout
            _selectedCellClusterInfoStatusBar->addWidget(descriptionLabel);


            std::vector<ClusterOrderContainer> orderedClustersSet;

            for (const auto& cluster : clusterValues) {
                ClusterOrderContainer temp{
                    static_cast<int>(cluster.getIndices().size()),
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
                auto clusterLabel = new ClickableLabel(); // Create the label without text
                QString labelText = QString("%1: %2").arg(clustersFromSet.name).arg(clustersFromSet.count);
                clusterLabel->setText(labelText); // Set the text on the label

                QColor textColor = clustersFromSet.color.lightness() > 127 ? Qt::black : Qt::white;
                clusterLabel->setStyleSheet(QString("ClickableLabel { color: %1; background-color: %2; padding: 2px; border: 0.5px solid %3; }")
                    .arg(textColor.name()).arg(clustersFromSet.color.name(QColor::HexArgb)).arg(textColor.name()));
                connect(clusterLabel, &ClickableLabel::clicked, this, [this, clusterLabel]() {


                    int current = _clusterCountSortingType.getCurrentIndex();
                    int newIndex;
                    if (current == 0)
                    {
                        newIndex = 1;
                    }
                    else if (current == 1)
                    {

                        if (!_customOrderClustersFromHierarchy.empty())
                        {
                            newIndex = 2;
                        }
                        else
                        {
                            newIndex = 0;
                        }
                    }
                    else
                    {
                        newIndex = 0;
                    }
                    _clusterCountSortingType.setCurrentIndex(newIndex);
                    });

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
}


void SettingsAction::setModifiedTriggeredData(QVariant geneListTable)
{
    if (!geneListTable.isNull())
    {
        ////startCodeTimer("Part15");
        //_filteredGeneNamesVariant.setVariant(geneListTable);
        _listModel.setVariant(geneListTable);
        ////stopCodeTimer("Part15");

    }
    else
    {
        qDebug() << "QVariant empty";
    }
}

void createTreeInitial(QJsonObject& node, const std::map<QString, InitialStatistics>& utilityMap) {
    // Check if the "name" key exists in the current node
    if (node.contains("name")) {
        QString nodeName = node["name"].toString();
        auto it = utilityMap.find(nodeName);

        if (it != utilityMap.end()) {
            node["mean"] = std::round(it->second.meanVal * 100.0) / 100.0; // Round to 2 decimal places
            node["differential"] = std::round(it->second.differentialVal * 100.0) / 100.0; // Round to 2 decimal places
            node["abundance"] = it->second.abundanceVal;
            node["rank"] = it->second.rankVal;
        }
    }

    // If the node has "children", recursively update them as well
    if (node.contains("children")) {
        QJsonArray children = node["children"].toArray();
        for (int i = 0; i < children.size(); ++i) {
            QJsonObject child = children[i].toObject();
            createTreeInitial(child, utilityMap); // Recursive call
            children[i] = child; // Update the modified object back into the array
        }
        node["children"] = children; // Update the modified array back into the parent JSON object
    }
}


void SettingsAction::precomputeTreesFromHierarchy()
{
    if (_mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker.isChecked())
    {
        return;
    }
    _precomputedTreesFromTheHierarchy.clear();
    auto start = std::chrono::high_resolution_clock::now();
    qDebug() << "Computing precomputeTreesFromHierarchy";

    if (!_speciesNamesDataset.getCurrentDataset().isValid() || !_mainPointsDataset.getCurrentDataset().isValid() || !_topClusterNamesDataset.getCurrentDataset().isValid() || !_middleClusterNamesDataset.getCurrentDataset().isValid() || !_bottomClusterNamesDataset.getCurrentDataset().isValid() || !_referenceTreeDataset.getCurrentDataset().isValid()) {
        qDebug() << "Datasets are not valid";
        return;
    }
    auto speciesNamesDataset = mv::data().getDataset<Clusters>(_speciesNamesDataset.getCurrentDataset().getDatasetId());
    auto mainPointsDataset = mv::data().getDataset<Points>(_mainPointsDataset.getCurrentDataset().getDatasetId());
    auto topClusterNamesDataset = mv::data().getDataset<Clusters>(_topClusterNamesDataset.getCurrentDataset().getDatasetId());
    auto middleClusterNamesDataset = mv::data().getDataset<Clusters>(_middleClusterNamesDataset.getCurrentDataset().getDatasetId());
    auto bottomClusterNamesDataset = mv::data().getDataset<Clusters>(_bottomClusterNamesDataset.getCurrentDataset().getDatasetId());

    auto referenceTreeDataset = mv::data().getDataset<CrossSpeciesComparisonTree>(_referenceTreeDataset.getCurrentDataset().getDatasetId());
    QJsonObject speciesDataJson = referenceTreeDataset->getTreeData();
    QStringList speciesNamesVerify = referenceTreeDataset->getTreeLeafNames();
    if (speciesDataJson.isEmpty() || speciesNamesVerify.isEmpty())
    {
        return;
    }


    if (speciesNamesDataset.isValid() && mainPointsDataset.isValid() && topClusterNamesDataset.isValid() && middleClusterNamesDataset.isValid() && bottomClusterNamesDataset.isValid())
    {
        auto speciesClusters = speciesNamesDataset->getClusters();

        //check if speciesNamesVerify and speciesClusters contain the same strings maybe in different order but same number and same value
        if (speciesNamesVerify.size() != speciesClusters.size())
        {
            return;
        }
        for (int i = 0; i < speciesNamesVerify.size(); i++)
        {
            if (speciesNamesVerify[i] != speciesClusters[i].getName())
            {
                return;
            }
        }


        auto mainPointDimensionNames = mainPointsDataset->getDimensionNames();
        auto mainPointsNumOfIndices = mainPointsDataset->getNumPoints();
        auto mainPointsNumOfDims = mainPointsDataset->getNumDimensions();

        QVector<Cluster> topClusters = topClusterNamesDataset->getClusters();
        QVector<Cluster>  middleClusters = middleClusterNamesDataset->getClusters();
        QVector<Cluster>  bottomClusters = bottomClusterNamesDataset->getClusters();

        if (!mainPointDimensionNames.empty()) {
            std::map<QString, QVector<Cluster>> combinedClusters = {
                {"top", topClusters},
                {"middle", middleClusters},
                {"bottom", bottomClusters}
            };

            QMutex mutex; // Mutex for thread safety

            QtConcurrent::blockingMap(combinedClusters, [&](const auto& pair) {
                const auto& hierarchyType = pair.first;
                const auto& clusters = pair.second;

                for (const auto& cluster : clusters) {
                    const auto& clusterName = cluster.getName();
                    auto clusterIndices = cluster.getIndices();
                    std::sort(clusterIndices.begin(), clusterIndices.end());
                    std::map<QString, std::map<QString, Stats>> topSpeciesToGeneExpressionMap;

                    QtConcurrent::blockingMap(speciesClusters, [&](const auto& species) {
                        const auto& speciesName = species.getName();
                        auto speciesIndices = species.getIndices();
                        std::sort(speciesIndices.begin(), speciesIndices.end());

                        std::vector<int> commonPointsIndices;
                        std::set_intersection(speciesIndices.begin(), speciesIndices.end(), clusterIndices.begin(), clusterIndices.end(), std::back_inserter(commonPointsIndices));

                        if (commonPointsIndices.empty()) {
                            return;
                        }

                        std::vector<float> resultContainerShort(commonPointsIndices.size());
                        std::vector<int> geneIndexContainer(1);

                        QtConcurrent::blockingMap(mainPointDimensionNames, [&](const QString& geneName) {
                            auto it = std::find(mainPointDimensionNames.begin(), mainPointDimensionNames.end(), geneName);
                            int geneIndex = (it != mainPointDimensionNames.end()) ? std::distance(mainPointDimensionNames.begin(), it) : -1;
                            geneIndexContainer[0] = geneIndex;

                            const auto& nonSelectionDetails = _clusterGeneMeanExpressionMap[speciesName][geneName]["allCells"];
                            int allCellCounts = nonSelectionDetails.first;
                            float allCellMean = nonSelectionDetails.second;

                            mainPointsDataset->populateDataForDimensions(resultContainerShort, geneIndexContainer, commonPointsIndices);

                            StatisticsSingle calculateStatisticsShort = calculateStatistics(resultContainerShort);

                            float allCellTotal = allCellMean * allCellCounts;
                            int nonSelectedCells = allCellCounts - calculateStatisticsShort.countVal;
                            float nonSelectedMean = (nonSelectedCells > 0) ? (allCellTotal - calculateStatisticsShort.meanVal * calculateStatisticsShort.countVal) / nonSelectedCells : 0.0f;

                            StatisticsSingle calculateStatisticsNot = { nonSelectedMean, nonSelectedCells };
                            int topHierarchyCountValue = (_clusterSpeciesFrequencyMap.find(speciesName) != _clusterSpeciesFrequencyMap.end()) ? _clusterSpeciesFrequencyMap[speciesName]["topCells"] : 0;
                            float topHierarchyFrequencyValue = (topHierarchyCountValue != 0) ? static_cast<float>(calculateStatisticsShort.countVal) / topHierarchyCountValue : 0.0f;

                            QMutexLocker locker(&mutex); // Lock the mutex for thread safety
                            topSpeciesToGeneExpressionMap[speciesName][geneName] = combineStatisticsSingle(calculateStatisticsShort, calculateStatisticsNot, topHierarchyCountValue);
                        });
                    });

                    enum class SelectionOption {
                        AbsoluteTopN,
                        PositiveTopN,
                        NegativeTopN
                    };

                    auto optionValue = "Positive";
                    SelectionOption option = SelectionOption::AbsoluteTopN;
                    if (optionValue == "Positive") {
                        option = SelectionOption::PositiveTopN;
                    } else if (optionValue == "Negative") {
                        option = SelectionOption::NegativeTopN;
                    }

                    std::map<QString, std::vector<std::pair<QString, int>>> rankingMap;

                    for (const auto& [speciesName, geneMap] : topSpeciesToGeneExpressionMap) {
                        std::vector<std::pair<QString, float>> geneExpressionVec;
                        geneExpressionVec.reserve(geneMap.size());
                        for (const auto& [geneName, stats] : geneMap) {
                            float differenceMeanValue = stats.meanSelected - stats.meanNonSelected;
                            geneExpressionVec.emplace_back(geneName, differenceMeanValue);
                        }

                        if (option == SelectionOption::AbsoluteTopN) {
                            std::sort(geneExpressionVec.begin(), geneExpressionVec.end(), [](const auto& a, const auto& b) {
                                return std::abs(a.second) > std::abs(b.second);
                            });
                        } else {
                            std::sort(geneExpressionVec.begin(), geneExpressionVec.end(), [](const auto& a, const auto& b) {
                                return a.second > b.second;
                            });
                            if (option == SelectionOption::NegativeTopN) {
                                std::reverse(geneExpressionVec.begin(), geneExpressionVec.end());
                            }
                        }

                        for (int i = 0; i < geneExpressionVec.size(); ++i) {
                            int rank = (option == SelectionOption::NegativeTopN) ? geneExpressionVec.size() - i : i + 1;
                            QMutexLocker locker(&mutex); // Lock the mutex for thread safety
                            rankingMap[geneExpressionVec[i].first].emplace_back(speciesName, rank);
                        }
                    }

                    for (auto& [geneName, speciesRankVec] : rankingMap) {
                        std::map<QString, InitialStatistics> utilityMap;
                        for (const auto& [speciesName, rank] : speciesRankVec) {
                            InitialStatistics tempStats;
                            tempStats.rankVal = rank;
                            tempStats.meanVal = topSpeciesToGeneExpressionMap[speciesName][geneName].meanSelected;
                            tempStats.differentialVal = topSpeciesToGeneExpressionMap[speciesName][geneName].meanSelected - topSpeciesToGeneExpressionMap[speciesName][geneName].meanNonSelected;
                            tempStats.abundanceVal = (topSpeciesToGeneExpressionMap[speciesName][geneName].abundanceCountTop != 0) ? topSpeciesToGeneExpressionMap[speciesName][geneName].meanSelected / topSpeciesToGeneExpressionMap[speciesName][geneName].abundanceCountTop : 0.0f;
                            utilityMap[speciesName] = tempStats;
                        }

                        QMutexLocker locker(&mutex); // Lock the mutex for thread safety
                        createTreeInitial(speciesDataJson, utilityMap);
                        _precomputedTreesFromTheHierarchy[hierarchyType][clusterName][geneName] = speciesDataJson;
                    }
                }
            });
        }



    }
    else
    {
        qDebug() << "Datasets are not valid";
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    qDebug() << "Time taken for precomputeTreesFromHierarchy : " + QString::number(duration / 1000.0) + " s";

}


void SettingsAction::computeGeneMeanExpressionMap()
{
    if (_mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker.isChecked())
    {
        return;
    }


    _clusterGeneMeanExpressionMap.clear();
    auto start = std::chrono::high_resolution_clock::now();
    qDebug() << "Computing gene mean expression map";

    _clusterGeneMeanExpressionMap.clear();
    if (_speciesNamesDataset.getCurrentDataset().isValid() && _mainPointsDataset.getCurrentDataset().isValid()) {
        auto speciesClusterDatasetFull = mv::data().getDataset<Clusters>(_speciesNamesDataset.getCurrentDataset().getDatasetId());
        auto mainPointDatasetFull = mv::data().getDataset<Points>(_mainPointsDataset.getCurrentDataset().getDatasetId());
        if (speciesClusterDatasetFull.isValid() && mainPointDatasetFull.isValid()) {
            auto speciesclusters = speciesClusterDatasetFull->getClusters();
            auto mainPointDimensionNames = mainPointDatasetFull->getDimensionNames();

            QtConcurrent::blockingMap(speciesclusters, [&](const auto& species) {
                auto speciesIndices = species.getIndices();
                auto speciesName = species.getName();
                for (int i = 0; i < mainPointDimensionNames.size(); i++) {
                    auto& geneName = mainPointDimensionNames[i];
                    auto geneIndex = { i };
                    std::vector<float> resultContainerFull(speciesIndices.size());
                    mainPointDatasetFull->populateDataForDimensions(resultContainerFull, geneIndex, speciesIndices);
                    float fullMean = calculateMean(resultContainerFull);
                    _clusterGeneMeanExpressionMap[speciesName][geneName]["allCells"] = std::make_pair(speciesIndices.size(), fullMean);
                }
                });

            _meanMapComputed = true;
        }
    }


    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    qDebug() << "Time taken for computeGeneMeanExpressionMap : " + QString::number(duration / 1000.0) + " s";
}


void SettingsAction::computeFrequencyMapForHierarchyItemsChange(QString hierarchyType)
{
    if (_mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker.isChecked() || hierarchyType.isEmpty())
    {
        return;
    }

    auto startTimer = std::chrono::high_resolution_clock::now();
    qDebug() << "computeFrequencyMapForHierarchyItemsChange Start for " + hierarchyType;

    if (!_speciesNamesDataset.getCurrentDataset().isValid() || !_mainPointsDataset.getCurrentDataset().isValid()) {
        return;
    }

    auto speciesClusterDatasetFull = mv::data().getDataset<Clusters>(_speciesNamesDataset.getCurrentDataset().getDatasetId());
    auto mainPointDatasetFull = mv::data().getDataset<Points>(_mainPointsDataset.getCurrentDataset().getDatasetId());
    auto numOfPoints = mainPointDatasetFull->getNumPoints();
    std::vector<bool> clusterNames(numOfPoints, true);
    QStringList inclusionList;
    mv::Dataset<Clusters> clusterDataset;

    if (hierarchyType == "top" && _topClusterNamesDataset.getCurrentDataset().isValid())
    {
        inclusionList = _topHierarchyClusterNamesFrequencyInclusionList.getSelectedOptions();
        clusterDataset = mv::data().getDataset<Clusters>(_topClusterNamesDataset.getCurrentDataset().getDatasetId());
    }
    else if (hierarchyType == "middle" && _middleClusterNamesDataset.getCurrentDataset().isValid())
    {
        inclusionList = _middleHierarchyClusterNamesFrequencyInclusionList.getSelectedOptions();
        clusterDataset = mv::data().getDataset<Clusters>(_middleClusterNamesDataset.getCurrentDataset().getDatasetId());
    }
    else if (hierarchyType == "bottom" && _bottomClusterNamesDataset.getCurrentDataset().isValid())
    {
        inclusionList = _bottomHierarchyClusterNamesFrequencyInclusionList.getSelectedOptions();
        clusterDataset = mv::data().getDataset<Clusters>(_bottomClusterNamesDataset.getCurrentDataset().getDatasetId());
    }

    if (clusterDataset.isValid())
    {
        for (const auto& cluster : clusterDataset->getClusters())
        {
            if (!inclusionList.contains(cluster.getName()))
            {
                for (const auto& index : cluster.getIndices())
                {
                    clusterNames[index] = false;
                }
            }
        }
    }

    if (speciesClusterDatasetFull.isValid() && mainPointDatasetFull.isValid())
    {
        auto speciesclusters = speciesClusterDatasetFull->getClusters();
        for (const auto& species : speciesclusters) {
            auto speciesIndices = species.getIndices();
            auto speciesName = species.getName();
            int count = std::count_if(speciesIndices.begin(), speciesIndices.end(), [&clusterNames](int index) {
                return clusterNames[index];
                });

            if (hierarchyType == "top")
            {
                _clusterSpeciesFrequencyMap[speciesName]["topCells"] = count;
            }
            else if (hierarchyType == "middle")
            {
                _clusterSpeciesFrequencyMap[speciesName]["middleCells"] = count;
            }
            else if (hierarchyType == "bottom")
            {
                _clusterSpeciesFrequencyMap[speciesName]["bottomCells"] = count;
            }
        }
    }

    auto endTimer = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTimer - startTimer).count();
    qDebug() << "Time taken for computeFrequencyMapForHierarchyItemsChange for " + hierarchyType + " : " + QString::number(duration / 1000.0) + " s";
}

void SettingsAction::computeGeneMeanExpressionMapForHierarchyItemsChangeExperimental(QString hierarchyType)
{
    if (_mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker.isChecked())
    {
        return;
    }
    auto startTimer = std::chrono::high_resolution_clock::now();
    qDebug() << "computeGeneMeanExpressionMapForHierarchyItemsChange Experimental Start for " + hierarchyType;
    if (hierarchyType == "")
    {
        return;
    }


    if (_speciesNamesDataset.getCurrentDataset().isValid() && _mainPointsDataset.getCurrentDataset().isValid()) {

        auto speciesClusterDatasetFull = mv::data().getDataset<Clusters>(_speciesNamesDataset.getCurrentDataset().getDatasetId());
        auto mainPointDatasetFull = mv::data().getDataset<Points>(_mainPointsDataset.getCurrentDataset().getDatasetId());
        auto numOfPoints = mainPointDatasetFull->getNumPoints();
        std::vector<bool> clusterNames(numOfPoints, true);
        bool datasetCheck = false;
        QStringList inclusionList;
        if (hierarchyType == "top")
        {
            inclusionList = _topHierarchyClusterNamesFrequencyInclusionList.getSelectedOptions();
            if (_topClusterNamesDataset.getCurrentDataset().isValid())
            {
                datasetCheck = true;
            }
        }
        else if (hierarchyType == "middle")
        {
            inclusionList = _middleHierarchyClusterNamesFrequencyInclusionList.getSelectedOptions();
            if (_middleClusterNamesDataset.getCurrentDataset().isValid())
            {
                datasetCheck = true;
            }
        }
        else if (hierarchyType == "bottom")
        {
            inclusionList = _bottomHierarchyClusterNamesFrequencyInclusionList.getSelectedOptions();
            if (_bottomClusterNamesDataset.getCurrentDataset().isValid())
            {
                datasetCheck = true;
            }
        }

        if (datasetCheck)
        {
            mv::Dataset<Clusters> clusterDataset;

            if (hierarchyType == "top")
            {
                clusterDataset = mv::data().getDataset<Clusters>(_topClusterNamesDataset.getCurrentDataset().getDatasetId());
            }
            else if (hierarchyType == "middle")
            {
                clusterDataset = mv::data().getDataset<Clusters>(_middleClusterNamesDataset.getCurrentDataset().getDatasetId());
            }
            else if (hierarchyType == "bottom")
            {
                clusterDataset = mv::data().getDataset<Clusters>(_bottomClusterNamesDataset.getCurrentDataset().getDatasetId());
            }

            for (auto cluster : clusterDataset->getClusters())
            {
                auto clusterIndices = cluster.getIndices();
                auto clusterName = cluster.getName();
                if (!inclusionList.contains(clusterName))
                {
                    for (auto index : clusterIndices)
                    {
                        clusterNames[index] = false;
                    }
                }

            }
        }



        if (speciesClusterDatasetFull.isValid() && mainPointDatasetFull.isValid())
        {
            auto speciesclusters = speciesClusterDatasetFull->getClusters();
            auto mainPointDimensionNames = mainPointDatasetFull->getDimensionNames();
            for (auto species : speciesclusters) {
                auto speciesIndices = species.getIndices();
                auto speciesName = species.getName();
                std::vector<int> indices;

                // Loop through all species indices to determine if they are in respective clusters
                for (int i = 0; i < speciesIndices.size(); ++i) {
                    // Check if the current species index is present in the cluster and only include those that are true
                    if (std::find(clusterNames.begin(), clusterNames.end(), speciesIndices[i]) != clusterNames.end()) {
                        indices.push_back(i);
                    }

                }

                for (int i = 0; i < mainPointDimensionNames.size(); i++) {
                    auto& geneName = mainPointDimensionNames[i];
                    auto geneIndex = { i };



                    std::vector<float> resultContainer(indices.size());
                    mainPointDatasetFull->populateDataForDimensions(resultContainer, geneIndex, indices);
                    float topMean = calculateMean(resultContainer);

                    if (hierarchyType == "top")
                    {
                        _clusterGeneMeanExpressionMap[speciesName][geneName]["topCells"] = std::make_pair(indices.size(), topMean);
                    }
                    else if (hierarchyType == "middle")
                    {
                        _clusterGeneMeanExpressionMap[speciesName][geneName]["middleCells"] = std::make_pair(indices.size(), topMean);
                    }
                    else if (hierarchyType == "bottom")
                    {
                        _clusterGeneMeanExpressionMap[speciesName][geneName]["bottomCells"] = std::make_pair(indices.size(), topMean);
                    }
                }

            }


        }
    }
    auto endTimer = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTimer - startTimer).count();
    qDebug() << "Time taken for computeGeneMeanExpressionMapForHierarchyItemsChangeExperimental for " + hierarchyType + " : " + QString::number(duration / 1000.0) + " s";

}
void SettingsAction::computeGeneMeanExpressionMapExperimental()
{
    if (_mapForHierarchyItemsChangeMethodStopForProjectLoadBlocker.isChecked())
    {
        return;
    }
    auto start = std::chrono::high_resolution_clock::now();
    qDebug() << "Computing gene mean expression map";


    _clusterGeneMeanExpressionMap.clear();

    if (_speciesNamesDataset.getCurrentDataset().isValid() && _mainPointsDataset.getCurrentDataset().isValid()) {

        auto speciesClusterDatasetFull = mv::data().getDataset<Clusters>(_speciesNamesDataset.getCurrentDataset().getDatasetId());
        auto mainPointDatasetFull = mv::data().getDataset<Points>(_mainPointsDataset.getCurrentDataset().getDatasetId());
        auto numOfPoints = mainPointDatasetFull->getNumPoints();
        std::vector<bool> topClusterNames(numOfPoints, true);
        std::vector<bool> middleClusterNames(numOfPoints, true);
        std::vector<bool>  bottomClusterNames(numOfPoints, true);
        QStringList topInclusionList = _topHierarchyClusterNamesFrequencyInclusionList.getSelectedOptions();
        QStringList middleInclusionList = _middleHierarchyClusterNamesFrequencyInclusionList.getSelectedOptions();
        QStringList bottomInclusionList = _bottomHierarchyClusterNamesFrequencyInclusionList.getSelectedOptions();
        if (_topClusterNamesDataset.getCurrentDataset().isValid() && _middleClusterNamesDataset.getCurrentDataset().isValid() && _bottomClusterNamesDataset.getCurrentDataset().isValid())

        {
            auto topClusterDataset = mv::data().getDataset<Clusters>(_topClusterNamesDataset.getCurrentDataset().getDatasetId());
            auto middleClusterDataset = mv::data().getDataset<Clusters>(_middleClusterNamesDataset.getCurrentDataset().getDatasetId());
            auto bottomClusterDataset = mv::data().getDataset<Clusters>(_bottomClusterNamesDataset.getCurrentDataset().getDatasetId());

            auto processCluster = [&](const Clusters& dataset, std::vector<bool>& clusterNames) {
                for (const auto& cluster : dataset.getClusters()) {
                    auto clusterIndices = cluster.getIndices();
                    auto clusterName = cluster.getName();
                    if (!topInclusionList.contains(clusterName)) {
                        for (auto index : clusterIndices) {
                            if (index < clusterNames.size()) {
                                clusterNames[index] = false;
                            }
                        }
                    }
                }
                };

            QFuture<void> topFuture = QtConcurrent::run([&]() { processCluster(*topClusterDataset, topClusterNames); });
            QFuture<void> middleFuture = QtConcurrent::run([&]() { processCluster(*middleClusterDataset, middleClusterNames); });
            QFuture<void> bottomFuture = QtConcurrent::run([&]() { processCluster(*bottomClusterDataset, bottomClusterNames); });

            topFuture.waitForFinished();
            middleFuture.waitForFinished();
            bottomFuture.waitForFinished();
        }



        // Ensure that the types match
        QMutex mapMutex; // Mutex to protect shared access to _clusterGeneMeanExpressionMap

        if (speciesClusterDatasetFull.isValid() && mainPointDatasetFull.isValid()) {
            auto speciesclusters = speciesClusterDatasetFull->getClusters();
            auto mainPointDimensionNames = mainPointDatasetFull->getDimensionNames();

            // Parallel processing of species clusters
            QtConcurrent::blockingMap(speciesclusters, [&](const auto& species) {
                auto speciesIndices = species.getIndices();
                auto speciesName = species.getName();

                std::vector<uint32_t> topIndices;
                std::vector<uint32_t> middleIndices;
                std::vector<uint32_t> bottomIndices;

                // Determine cluster membership for the species
                for (uint32_t i = 0; i < speciesIndices.size(); ++i) {
                    if (std::binary_search(topClusterNames.begin(), topClusterNames.end(), speciesIndices[i])) {
                        topIndices.push_back(i);
                    }
                    if (std::binary_search(middleClusterNames.begin(), middleClusterNames.end(), speciesIndices[i])) {
                        middleIndices.push_back(i);
                    }
                    if (std::binary_search(bottomClusterNames.begin(), bottomClusterNames.end(), speciesIndices[i])) {
                        bottomIndices.push_back(i);
                    }
                }

                // Parallel processing of gene expressions within each species
                QtConcurrent::blockingMap(mainPointDimensionNames, [&](const auto& geneName) {
                    // Manually find the index of the geneName
                    int geneIndex = std::distance(mainPointDimensionNames.begin(),
                    std::find(mainPointDimensionNames.begin(), mainPointDimensionNames.end(), geneName));

                if (geneIndex == mainPointDimensionNames.size()) {
                    // Handle case where geneName is not found if necessary
                    return; // Skip processing if the index is invalid
                }

                auto processCells = [&](const std::vector<uint32_t>& indices, const QString& cellType) {
                    std::vector<float> resultContainer(indices.size());
                    mainPointDatasetFull->populateDataForDimensions(resultContainer, std::vector<int>{geneIndex}, indices);
                    float mean = calculateMean(resultContainer);
                    QMutexLocker locker(&mapMutex);
                    _clusterGeneMeanExpressionMap[speciesName][geneName][cellType] = std::make_pair(indices.size(), mean);
                    };

                processCells(speciesIndices, "allCells");
                processCells(topIndices, "topCells");
                processCells(middleIndices, "middleCells");
                processCells(bottomIndices, "bottomCells");
                    });
                });

            _meanMapComputed = true;
        }



    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    qDebug() << "\n\n++++++++++++++++++Time taken for computeGeneMeanExpressionMap : " + QString::number(duration / 1000.0) + " s";

}
void SettingsAction::findTopNGenesPerCluster() {

    int n = _topNGenesFilter.getValue();

    if (_clusterNameToGeneNameToExpressionValue.empty() || n <= 0) {
        return;
    }

    // startCodeTimer("findTopNGenesPerCluster");

    enum class SelectionOption {
        AbsoluteTopN,
        PositiveTopN,
        NegativeTopN
    };
    auto optionValue = _typeofTopNGenes.getCurrentText();
    SelectionOption option = SelectionOption::AbsoluteTopN;
    if (optionValue == "Positive") {
        option = SelectionOption::PositiveTopN;
    }
    else if (optionValue == "Negative") {
        option = SelectionOption::NegativeTopN;
    }

    _uniqueReturnGeneList.clear();
    std::map<QString, std::vector<QString>> geneAppearanceCounter;
    std::map<QString, std::vector<std::pair<QString, int>>> rankingMap;

    for (const auto& outerPair : _clusterNameToGeneNameToExpressionValue) {
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
                    _uniqueReturnGeneList.insert(geneExpressionVec[i].first);
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
                    _uniqueReturnGeneList.insert(geneExpressionVec[i].first);
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
                    _uniqueReturnGeneList.insert(geneExpressionVec[i].first);
                    if (outerPair.second.find(geneExpressionVec[i].first)->second.meanSelected > 0) {
                        geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                    }

                    //geneAppearanceCounter[geneExpressionVec[i].first].push_back(speciesName);
                }
                rankingMap[geneExpressionVec[i].first].emplace_back(speciesName, geneExpressionVec.size() - i); // Corrected rank calculation
            }
            break;
        }
        }
    }

    //iterate std::map<QString, std::vector<std::pair<QString, int>>> rankingMap;
    // Iterating over the map
    if (_performGeneTableTsneAction.isChecked()) {
        std::vector<float> rankOrder;

        std::unordered_map<QString, std::unordered_map<QString, float>> geneSimilarityMap;
        std::unordered_set<QString> speciesSet;
        std::unordered_set<QString> geneSet;

        for (const auto& item : rankingMap) {
            const QString& gene = item.first;
            std::unordered_map<QString, float> rankCounter;
            for (const auto& pair : item.second) {
                const QString& species = pair.first;
                const float rank = (pair.second <= n) ? 1.0f : 0.0f; // Use float directly
                rankCounter[species] = rank;
                speciesSet.insert(species);
            }
            geneSimilarityMap[gene] = std::move(rankCounter);
            geneSet.insert(gene);
        }
        _geneOrder.clear();
        std::vector<QString> speciesOrder(speciesSet.begin(), speciesSet.end());
        _geneOrder = std::vector<QString>(geneSet.begin(), geneSet.end());
        rankOrder.resize(_geneOrder.size() * speciesOrder.size(), 0.0f); // Initialize with 0.0f for clarity

        std::unordered_map<QString, int> speciesIndexMap;
        for (int i = 0; i < speciesOrder.size(); ++i) {
            speciesIndexMap[speciesOrder[i]] = i;
        }

        for (int geneIndex = 0; geneIndex < _geneOrder.size(); ++geneIndex) {
            const QString& gene = _geneOrder[geneIndex];
            const auto& speciesRanks = geneSimilarityMap[gene];
            for (const auto& speciesRank : speciesRanks) {
                const QString& species = speciesRank.first;
                const float rank = speciesRank.second; // Already a float, no need to cast
                int speciesIndex = speciesIndexMap[species];
                rankOrder[geneIndex * speciesOrder.size() + speciesIndex] = rank;
            }
        }

    QString pointDataId = _geneSimilarityPoints->getId();
    int pointDimSize = static_cast<int>(speciesOrder.size());
    int pointIndicesSize = static_cast<int>(_geneOrder.size());

    if (_selectedPointsTSNEDatasetForGeneTable.isValid())
    {
        auto runningAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDatasetForGeneTable->findChildByPath("TSNE/TsneComputationAction/Running"));

        if (runningAction)
        {

            if (runningAction->isChecked())
            {
                auto stopAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDatasetForGeneTable->findChildByPath("TSNE/TsneComputationAction/Stop"));
                if (stopAction)
                {
                    stopAction->trigger();
                    //std::this_thread::sleep_for(std::chrono::seconds(5));
                }
            }

        }
        mv::data().removeDataset(_selectedPointsTSNEDatasetForGeneTable);
    }

    populatePointData(pointDataId, rankOrder, pointIndicesSize, pointDimSize, speciesOrder);

    mv::plugin::AnalysisPlugin* analysisPlugin;
    auto scatterplotModificationsGeneSimilarity = [this]() {
        if (_selectedPointsTSNEDatasetForGeneTable.isValid()) {
            auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
            mv::gui::DatasetPickerAction* colorDatasetPickerAction;
            mv::gui::DatasetPickerAction* pointDatasetPickerAction;
            mv::gui::ViewPluginSamplerAction* samplerActionAction;
            if (scatterplotViewFactory) {
                for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                    if (plugin->getGuiName() == "Scatterplot Cell Selection Overview") {
                        pointDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Position"));
                        if (pointDatasetPickerAction) {
                            pointDatasetPickerAction->setCurrentText("");

                            pointDatasetPickerAction->setCurrentDataset(_selectedPointsTSNEDatasetForGeneTable);

                            colorDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Color"));
                            if (colorDatasetPickerAction)
                            {
                                colorDatasetPickerAction->setCurrentText("");

                                if (_geneSimilarityClusterColoring.isValid())
                                {
                                    colorDatasetPickerAction->setCurrentDataset(_geneSimilarityClusterColoring);
                                }

                            }
                            
                            samplerActionAction = plugin->findChildByPath<mv::gui::ViewPluginSamplerAction>("Sampler");

                            if (samplerActionAction)
                            {
                                samplerActionAction->setTooltipGeneratorFunction([this](const ViewPluginSamplerAction::SampleContext& toolTipContext) -> QString {
                                    QString clusterDatasetId = _speciesNamesDataset.getCurrentDataset().getDatasetId();
                                    return generateTooltip(toolTipContext, clusterDatasetId,true, "GlobalPointIndices");
                                    });
                            }
                        }
                    }
                }
            }
        }

        };



    {
        //startCodeTimer("Part10");
        analysisPlugin = mv::plugins().requestPlugin<AnalysisPlugin>("tSNE Analysis", { _geneSimilarityPoints });
        if (!analysisPlugin) {
            qDebug() << "Could not find create TSNE Analysis";
            return;
        }
        _selectedPointsTSNEDatasetForGeneTable = analysisPlugin->getOutputDataset();
        int groupID2 = 10 * 3;
        _selectedPointsTSNEDatasetForGeneTable->setGroupIndex(groupID2);
        if (_selectedPointsTSNEDatasetForGeneTable.isValid())
        {
            //_selectedPointsTSNEDatasetForGeneTable->printChildren();
            bool skip = false;
            int perplexity = std::min(static_cast<int>(_geneOrder.size()), _tsnePerplexity.getValue());
            if (perplexity < 5)
            {
                qDebug() << "Perplexity is less than 5";
                skip = true;
                //_startComputationTriggerAction.setDisabled(false);
            }
            if (!skip)
            {
            if (perplexity != _tsnePerplexity.getValue())
            {
                _tsnePerplexity.setValue(perplexity);
            }

            auto perplexityAction = dynamic_cast<IntegralAction*>(_selectedPointsTSNEDatasetForGeneTable->findChildByPath("TSNE/Perplexity"));
            if (perplexityAction)
            {
                //qDebug() << "Perplexity: Found";
                perplexityAction->setValue(perplexity);
            }
            else
            {
                qDebug() << "Perplexity: Not Found";
            }

            QString knnAlgorithmValue = _performGeneTableTsneKnn.getCurrentText();
            QString distanceMetricValue = _performGeneTableTsneDistance.getCurrentText();
            if (knnAlgorithmValue != "")
            {
                auto knnAction = dynamic_cast<OptionAction*>(_selectedPointsTSNEDatasetForGeneTable->findChildByPath("TSNE/kNN Algorithm"));
                if (knnAction)
                {
                    //qDebug() << "Knn: Found";
                    try {
                        knnAction->setCurrentText(knnAlgorithmValue);
                    }
                    catch (const std::exception& e) {
                        qDebug() << "An exception occurred in setting knn value: " << e.what();
                    }  
                }
                else
                {
                    qDebug() << "Knn: Not Found";
                }
            }
            if (distanceMetricValue != "")
            {
                auto distanceAction = dynamic_cast<OptionAction*>(_selectedPointsTSNEDatasetForGeneTable->findChildByPath("TSNE/Distance metric"));
                if (distanceAction)
                {
                    //qDebug() << "Distance: Found";
                    try {
                        distanceAction->setCurrentText(distanceMetricValue);
                    }
                    catch (const std::exception& e) {
                        qDebug() << "An exception occurred in setting distance value: " << e.what();
                    }
                }
                else
                {
                    qDebug() << "Distance: Not Found";
                }
            }

            scatterplotModificationsGeneSimilarity();

            auto startAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDatasetForGeneTable->findChildByPath("TSNE/TsneComputationAction/Start"));
            if (startAction) {

                startAction->trigger();

                analysisPlugin->getOutputDataset()->setSelectionIndices({});
            }
        }
        }
        //stopCodeTimer("Part10");
    }

    std::vector<int> selectedIndices;
    std::vector<int> nonselectedIndices;
    selectedIndices.reserve(_geneOrder.size()); // Pre-allocate memory
    nonselectedIndices.reserve(_geneOrder.size()); // Pre-allocate memory
    for (int i = 0; i < _geneOrder.size(); i++)
    {
        if (_uniqueReturnGeneList.find(_geneOrder[i]) != _uniqueReturnGeneList.end())
        {
            selectedIndices.push_back(i);
        }
        else
        {
            nonselectedIndices.push_back(i);
        }
    }
    QString clusterDataId = _geneSimilarityClusterColoring->getId();
    QColor selectedColor = QColor("#00A2ED");
    QColor nonSelectedColor = QColor("#ff5d12");
    std::map<QString, std::pair<QColor, std::vector<int>>> selectedClusterMap;
    selectedClusterMap["TopNSelectedGenes"] = { selectedColor, selectedIndices };
    selectedClusterMap["NonTopNGenes"] = { nonSelectedColor, nonselectedIndices };

        populateClusterData(clusterDataId, selectedClusterMap);
}

    //stopCodeTimer("findTopNGenesPerCluster");
    QVariant returnedmodel= createModelFromData(_clusterNameToGeneNameToExpressionValue, geneAppearanceCounter, rankingMap, n);

    setModifiedTriggeredData(returnedmodel);

    //return returnedmodel;
}
void SettingsAction::removeDatasets(int groupId)
{
    auto allDatasets = mv::data().getAllDatasets();

    Datasets datasetsFilteredAndSorted;

    std::copy_if(allDatasets.begin(), allDatasets.end(), std::back_inserter(datasetsFilteredAndSorted), [groupId](Dataset<DatasetImpl> dataset) {
        return dataset->getGroupIndex() == groupId;
    });

    std::sort(datasetsFilteredAndSorted.begin(), datasetsFilteredAndSorted.end(), [](Dataset<DatasetImpl> lhs, Dataset<DatasetImpl> rhs) -> bool {
        return rhs->getDataHierarchyItem().getDepth() < lhs->getDataHierarchyItem().getDepth();
    });

    //std::reverse(datasetsFilteredAndSorted.begin(), datasetsFilteredAndSorted.end());

    for (auto dataset : datasetsFilteredAndSorted)
    {
        //qDebug() << dataset->getGuiName() << dataset->getId() << dataset->getDataHierarchyItem().getDepth();

        if (dataset.isValid())
            mv::data().removeDataset(dataset);
    }

    //mv::data().removeDataset(_selectedPointsDataset);
    //mv::data().removeDataset(_selectedPointsEmbeddingDataset);
}
QVariant SettingsAction::createModelFromData(const std::map<QString, std::map<QString, Stats>>& map,   const std::map<QString, std::vector<QString>>& geneCounter, const std::map<QString, std::vector<std::pair<QString, int>>>& rankingMap, const int& n) {

    if (map.empty() || _totalGeneList.empty()) {
        return QVariant();
    }
    //startCodeTimer("createModelFromData");
    QStandardItemModel* model = new QStandardItemModel();
    _initColumnNames = { "ID", "Species \nAppearance", "Gene Appearance Species Names", "Statistics" };
    model->setHorizontalHeaderLabels(_initColumnNames);

    QStringList headers = _initColumnNames;
    _hiddenShowncolumns.setOptions(headers);
    _hiddenShowncolumns.setSelectedOptions({ headers[0], headers[1] });

    for (const auto& gene : _totalGeneList) {
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

        row.push_back(new QStandardItem(gene)); // ID(string) should sort by string

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
            count = static_cast<int>(speciesDetails.size());
            QStringList speciesNames;
            for (const auto& speciesDetail : speciesDetails) {
                speciesNames << speciesDetail;
            }
            speciesGeneAppearancesComb = speciesNames.join(";");
        }
        auto* countItem = new QStandardItem(); // Gene Appearances (int) should sort by int
        countItem->setData(count, Qt::DisplayRole);
        countItem->setData(count, Qt::UserRole); // Use Qt::UserRole or another custom role for sorting by integer
        row.push_back(countItem);

        //row.push_back(new QStandardItem(QString::number(count))); // Gene Appearances (int) should sort by int
        row.push_back(new QStandardItem(speciesGeneAppearancesComb)); // Gene Appearance Species Names (string) should sort by string

        QString formattedStatistics;
        for (const auto& [species, stats] : statisticsValuesForSpeciesMap) {
            formattedStatistics += QString("Species: %1, Rank: %2, AbundanceTop: %3, MeanSelected: %4, CountSelected: %5, MeanNotSelected: %6, CountNotSelected: %7;\n")//, MeanAll: %7, CountAll: %8
                .arg(species)
                .arg(rankcounter[species])
                .arg(stats.abundanceCountTop)
                .arg(stats.meanSelected, 0, 'f', 2)
                .arg(stats.countSelected)
                .arg(stats.meanNonSelected, 0, 'f', 2)
                .arg(stats.countNonSelected)
                //.arg(stats.meanAll, 0, 'f', 2)
                //.arg(stats.countAll)
                ;
        }
        row.push_back(new QStandardItem(formattedStatistics)); // Statistics (string) should sort by string
        model->appendRow(row);
    }

    //stopCodeTimer("createModelFromData");

    return QVariant::fromValue(model);

}
void SettingsAction::createClusterPositionMap()
{
    _clusterPositionMap;
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
    (void)QtConcurrent::run([this, datasetId, pointVector, numPoints, numDimensions, dimensionNames]() {
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
    _toggleScatterplotSelection.setDisabled(false);
    _usePreComputedTSNE.setDisabled(false);
    _tsnePerplexity.setDisabled(false);
    _performGeneTableTsnePerplexity.setDisabled(false);
    _performGeneTableTsneKnn.setDisabled(false);
    _performGeneTableTsneDistance.setDisabled(false);
    _performGeneTableTsneTrigger.setDisabled(false);
    _referenceTreeDataset.setDisabled(false);
    _mainPointsDataset.setDisabled(false);
    _embeddingDataset.setDisabled(false);
    _speciesNamesDataset.setDisabled(false);
    _bottomClusterNamesDataset.setDisabled(false);
    _middleClusterNamesDataset.setDisabled(false);
    _topClusterNamesDataset.setDisabled(false);
    _speciesExplorerInMap.setDisabled(false);
    _topHierarchyClusterNamesFrequencyInclusionList.setDisabled(false);
    _middleHierarchyClusterNamesFrequencyInclusionList.setDisabled(false);
    _bottomHierarchyClusterNamesFrequencyInclusionList.setDisabled(false);
    _speciesExplorerInMapTrigger.setDisabled(false);
    _revertRowSelectionChangesToInitial.setDisabled(false);
    _scatterplotEmbeddingPointsUMAPOption.setDisabled(false);
    _selectedSpeciesVals.setDisabled(false);
    _clusterOrderHierarchy.setDisabled(false);
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
    _toggleScatterplotSelection.setChecked(true);
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
    _speciesExplorerInMapTrigger.setDisabled(true);
    _usePreComputedTSNE.setDisabled(true);
    _applyLogTransformation.setDisabled(true);
    _speciesExplorerInMap.setDisabled(true);
    _revertRowSelectionChangesToInitial.setDisabled(true);
    _toggleScatterplotSelection.setDisabled(true);
    _tsnePerplexity.setDisabled(true);
    _performGeneTableTsnePerplexity.setDisabled(true);
    _performGeneTableTsneKnn.setDisabled(true);
    _performGeneTableTsneDistance.setDisabled(true);
    _performGeneTableTsneTrigger.setDisabled(true);
    _referenceTreeDataset.setDisabled(true);
    _mainPointsDataset.setDisabled(true);
    _embeddingDataset.setDisabled(true);
    _speciesNamesDataset.setDisabled(true);
    _bottomClusterNamesDataset.setDisabled(true);
    _middleClusterNamesDataset.setDisabled(true);
    _topClusterNamesDataset.setDisabled(true);
    _scatterplotEmbeddingPointsUMAPOption.setDisabled(true);
    _topHierarchyClusterNamesFrequencyInclusionList.setDisabled(true);
    _middleHierarchyClusterNamesFrequencyInclusionList.setDisabled(true);
    _bottomHierarchyClusterNamesFrequencyInclusionList.setDisabled(true);
    _selectedSpeciesVals.setDisabled(true);
    _clusterOrderHierarchy.setDisabled(true);
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
    _revertRowSelectionChangesToInitial.setDisabled(false);
    _speciesExplorerInMapTrigger.setDisabled(false);
    //if (!stringActionHasOptions)
    //{
    //    _revertRowSelectionChangesToInitial.setDisabled(true);
    //}
    //else
    //{
    //    if (!optionsActionHasOptions)
    //    {

    //        _revertRowSelectionChangesToInitial.setDisabled(false);
    //    }
    //    else
    //    {
    //        if (bothListsEqual)
    //        {

    //            _revertRowSelectionChangesToInitial.setDisabled(true);
    //        }
    //        else
    //        {

    //            _revertRowSelectionChangesToInitial.setDisabled(false);
    //        }
    //    }
    //}



    //if (!optionsActionHasOptions)
    //{
    //    _speciesExplorerInMapTrigger.setDisabled(true);

    //}
    //else 
    //{
    //    _speciesExplorerInMapTrigger.setDisabled(false);
    //}



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
        colorDataset->getClusters() = QVector<Cluster>();
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

QString SettingsAction::generateTooltip(const ViewPluginSamplerAction::SampleContext& toolTipContext, const QString& clusterDatasetId, bool showTooltip, QString indicesType) {
    // Extract and convert GlobalPointIndices and ColorDatasetID from toolTipContext
    auto raw_Global_Local_PointIndices = toolTipContext[indicesType].toList();

    // Convert the list of global point indices to a vector of integers
    std::vector<std::seed_seq::result_type> global_local_PointIndices;
    global_local_PointIndices.reserve(raw_Global_Local_PointIndices.size());
    for (const auto& global_local_PointIndex : raw_Global_Local_PointIndices) {
        global_local_PointIndices.push_back(global_local_PointIndex.toInt());
    }

    // If the global point indices list is empty, return an empty result
    if (global_local_PointIndices.empty()) {
        return {};
    }

    // If there is no cluster dataset ID, return a summary of total points
    if (clusterDatasetId.isEmpty()) {
        return QString("<table> \
<tr> \
<td><b>Total points: </b></td> \
<td>%1</td> \
</tr> \
</table>").arg(global_local_PointIndices.size());
    }

    // Retrieve the cluster dataset
    auto clusterFullDataset = mv::data().getDataset<Clusters>(clusterDatasetId);

    // If the dataset is invalid, return a summary of total points
    if (!clusterFullDataset.isValid()) {
        return QString("<table> \
<tr> \
<td><b>Total points: </b></td> \
<td>%1</td> \
</tr> \
</table>").arg(global_local_PointIndices.size());
    }

    // Get the clusters from the dataset
    auto clusterValuesData = clusterFullDataset->getClusters();

    // If the clusters data is empty, return a summary of total points
    if (clusterValuesData.isEmpty()) {
        return QString("<table> \
<tr> \
<td><b>Total points: </b></td> \
<td>%1</td> \
</tr> \
</table>").arg(global_local_PointIndices.size());
    }

    // Process each cluster and find intersections with global point indices
    std::map<QString, std::pair<int, QColor>> clusterCountMap;
    for (const auto& cluster : clusterValuesData) {
        QString clusterName = cluster.getName();
        QColor clusterColor = cluster.getColor();
        const auto& clusterIndices = cluster.getIndices();

        std::vector<std::seed_seq::result_type> sortedClusterIndices = clusterIndices;
        std::vector<std::seed_seq::result_type> sortedGlobalLocalPointIndices = global_local_PointIndices;

        std::sort(sortedClusterIndices.begin(), sortedClusterIndices.end());
        std::sort(sortedGlobalLocalPointIndices.begin(), sortedGlobalLocalPointIndices.end());

        std::vector<std::seed_seq::result_type> intersect;
        intersect.reserve(std::min(sortedClusterIndices.size(), sortedGlobalLocalPointIndices.size()));
        std::set_intersection(sortedClusterIndices.begin(), sortedClusterIndices.end(),
            sortedGlobalLocalPointIndices.begin(), sortedGlobalLocalPointIndices.end(),
            std::back_inserter(intersect));

        // If there is an intersection, store the result in the map
        if (!intersect.empty()) {
            clusterCountMap[clusterName] = std::make_pair(intersect.size(), clusterColor);
        }
    }

    // If no clusters were found, return a summary of total points
    if (clusterCountMap.empty()) {
        return QString("<table> \
<tr> \
<td><b>Total points: </b></td> \
<td>%1</td> \
</tr> \
</table>").arg(global_local_PointIndices.size());
    }

    // Generate HTML output
    //QString html = "<html><head><style>"
    //    "table { border-collapse: collapse; width: 100%; font-size: 12px; }"
    //    "th, td { border: 1px solid black; padding: 4px; text-align: left; }"
    //    "th { background-color: #f2f2f2; }"
    //    "</style></head><body>";
    //html += "<table>";
    //html += "<tr><th>Cluster Name</th><th>Count</th></tr>";

    //// Populate the table with cluster data
    //for (const auto& entry : clusterCountMap) {
    //    QString clusterName = entry.first;
    //    int count = entry.second.first;
    //    QString colorHex = entry.second.second.name();

    //    html += "<tr>";
    //    html += "<td style='background-color:" + colorHex + ";'>" + clusterName + "</td>";
    //    html += "<td>" + QString::number(count) + "</td>";
    //    html += "</tr>";
    //}

    //html += "</table></body></html>";
    //return html;
    QString html = "<html><head><style>"
        "body { display: flex; flex-wrap: wrap; max-width: 800px; }" // Set max-width to control wrapping
        "div { display: inline-block; padding: 4px; margin: 2px; border: 1px solid black; font-size: 12px; }"
        "</style></head><body>";

    // Function to determine if a color is dark
    auto isDarkColor = [](const QColor& color) {
        int brightness = (color.red() * 299 + color.green() * 587 + color.blue() * 114) / 1000;
        return brightness < 128;
        };

    // Convert the map to a vector of pairs for sorting
    std::vector<std::pair<QString, std::pair<int, QColor>>> clusterVector(clusterCountMap.begin(), clusterCountMap.end());

    // Sort the vector by count in descending order
    std::sort(clusterVector.begin(), clusterVector.end(), [](const auto& a, const auto& b) {
        return a.second.first > b.second.first;
        });

    // Populate the divs with cluster data
    for (const auto& entry : clusterVector) {
        QString clusterName = entry.first;
        int count = entry.second.first;
        QString colorHex = entry.second.second.name();
        QColor color(entry.second.second);

        QString textColor = isDarkColor(color) ? "white" : "black";

        html += "<div style='background-color:" + colorHex + "; color:" + textColor + ";'>";
        html += clusterName + ":" + QString::number(count);
        html += "</div>";
    }

    html += "</body></html>";
    return html;




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
    _bottomClusterNamesDataset.fromParentVariantMap(variantMap);
    _middleClusterNamesDataset.fromParentVariantMap(variantMap);
    _topClusterNamesDataset.fromParentVariantMap(variantMap);
    _filteredGeneNamesVariant.fromParentVariantMap(variantMap);
    _topNGenesFilter.fromParentVariantMap(variantMap);
    _filteringEditTreeDataset.fromParentVariantMap(variantMap);
    _referenceTreeDataset.fromParentVariantMap(variantMap);
    _selectedRowIndex.fromParentVariantMap(variantMap);
    _performGeneTableTsneAction.fromParentVariantMap(variantMap);
    _tsnePerplexity.fromParentVariantMap(variantMap);
    _performGeneTableTsnePerplexity.fromParentVariantMap(variantMap);
    _performGeneTableTsneKnn.fromParentVariantMap(variantMap);
    _performGeneTableTsneDistance.fromParentVariantMap(variantMap);
    _performGeneTableTsneTrigger.fromParentVariantMap(variantMap);
    _hiddenShowncolumns.fromParentVariantMap(variantMap);
    _speciesExplorerInMap.fromParentVariantMap(variantMap);
    _topHierarchyClusterNamesFrequencyInclusionList.fromParentVariantMap(variantMap);
    _middleHierarchyClusterNamesFrequencyInclusionList.fromParentVariantMap(variantMap);
    _bottomHierarchyClusterNamesFrequencyInclusionList.fromParentVariantMap(variantMap);
    _scatterplotReembedColorOption.fromParentVariantMap(variantMap);
    _scatterplotEmbeddingPointsUMAPOption.fromParentVariantMap(variantMap);
    _selectedSpeciesVals.fromParentVariantMap(variantMap);
    _clusterOrderHierarchy.fromParentVariantMap(variantMap);
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
    _bottomClusterNamesDataset.insertIntoVariantMap(variantMap);
    _middleClusterNamesDataset.insertIntoVariantMap(variantMap);
    _topClusterNamesDataset.insertIntoVariantMap(variantMap);
    _filteredGeneNamesVariant.insertIntoVariantMap(variantMap);
    _topNGenesFilter.insertIntoVariantMap(variantMap);
    _filteringEditTreeDataset.insertIntoVariantMap(variantMap);
    _referenceTreeDataset.insertIntoVariantMap(variantMap);
    _selectedRowIndex.insertIntoVariantMap(variantMap);
    _performGeneTableTsneAction.insertIntoVariantMap(variantMap);
    _tsnePerplexity.insertIntoVariantMap(variantMap);
    _performGeneTableTsnePerplexity.insertIntoVariantMap(variantMap);
    _performGeneTableTsneDistance.insertIntoVariantMap(variantMap);
    _performGeneTableTsneKnn.insertIntoVariantMap(variantMap);
    _performGeneTableTsneTrigger.insertIntoVariantMap(variantMap);
    _hiddenShowncolumns.insertIntoVariantMap(variantMap);
    _speciesExplorerInMap.insertIntoVariantMap(variantMap);
    _topHierarchyClusterNamesFrequencyInclusionList.insertIntoVariantMap(variantMap);
    _middleHierarchyClusterNamesFrequencyInclusionList.insertIntoVariantMap(variantMap);
    _bottomHierarchyClusterNamesFrequencyInclusionList.insertIntoVariantMap(variantMap);
    _scatterplotReembedColorOption.insertIntoVariantMap(variantMap);
    _scatterplotEmbeddingPointsUMAPOption.insertIntoVariantMap(variantMap);
    _selectedSpeciesVals.insertIntoVariantMap(variantMap);
    _clusterOrderHierarchy.insertIntoVariantMap(variantMap);
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
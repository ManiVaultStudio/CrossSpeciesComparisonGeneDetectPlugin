#include "SettingsAction.h"
#include "CrossSpeciesComparisonGeneDetectPlugin.h"
#include<string>  
#include <QFileDialog>
#include <QPageLayout>
#include <QWebEngineView>
#include <CrossSpeciesComparisonTreeData.h>
#include <numeric>   // for std::reduce
#include <execution> // for std::execution::par
#include "lib/Distance/annoylib.h"
#include "lib/Distance/kissrandom.h"
#include <QtConcurrent>
#include "lib/JSONnlohmann/json.hpp"
#include "lib/Clustering/fastcluster.h"
#include "lib/NewickComparator/newick_comparator.h" //https://github.com/MaciejSurowiec/Maximum_agreement_subtree_problem
#include <sstream>
#include <stack>
#include <algorithm> // for std::find
#include <mutex>

#include <iostream>
#include <map>
#include <string>
#include <chrono>
#include <cmath> // Include for std::log

using namespace mv;
using namespace mv::gui;




std::map<std::string, std::chrono::high_resolution_clock::time_point> timers;

Statistics combineStatisticsSingle(const StatisticsSingle& selected, const StatisticsSingle& nonSelected) {
    Statistics combinedStats;
    combinedStats.meanSelected = selected.meanVal;
    combinedStats.medianSelected = selected.medianVal;
    combinedStats.modeSelected = selected.modeVal;
    combinedStats.rangeSelected = selected.rangeVal;
    combinedStats.countSelected = selected.countVal;

    combinedStats.meanNonSelected = nonSelected.meanVal; // Note the change here to meanSelected from nonSelected
    combinedStats.medianNonSelected = nonSelected.medianVal; // Same as above
    combinedStats.modeNonSelected = nonSelected.modeVal; // Same as above
    combinedStats.rangeNonSelected = nonSelected.rangeVal; // Same as above
    combinedStats.countNonSelected = nonSelected.countVal; // Same as above

    return combinedStats;
}
StatisticsSingle calculateStatistics(const std::vector<float>& data) {
    if (data.empty()) {
        return { std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(),
                 std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN() };
    }
    int count = data.size();
    // Calculate mean using std::reduce for better performance with parallel execution
    float sum = std::reduce(std::execution::par, data.begin(), data.end(), 0.0f);
    float mean = sum / data.size();

    // Sort data to find median and range
    std::vector<float> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());
    float median;
    size_t size = sortedData.size();
    if (size % 2 == 0) {
        median = (sortedData[size / 2 - 1] + sortedData[size / 2]) / 2;
    }
    else {
        median = sortedData[size / 2];
    }

    // Calculate mode
    std::unordered_map<float, int> frequency;
    for (float num : data) {
        frequency[num]++;
    }
    int maxCount = 0;
    float mode = std::numeric_limits<float>::quiet_NaN();
    for (const auto& pair : frequency) {
        if (pair.second > maxCount) {
            maxCount = pair.second;
            mode = pair.first;
        }
    }

    // Calculate range
    float range = sortedData.back() - sortedData.front();
    
    return { mean, median, mode, range, count };
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

float calculateMean(const std::vector<float>& v) {
    if (v.empty())
        return 0.0f;

    float sum = std::reduce(std::execution::par, v.begin(), v.end(), 0.0f);
    return sum / static_cast<float>(v.size());
}

float calculateMedian(const std::vector<float>& vec) {
    if (vec.empty()) return 0.0f;

    std::vector<float> v = vec; // Copy to avoid modifying the original vector
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + n, v.end()); // Partially sort to find the median

    if (v.size() % 2 == 1) {
        // For odd-sized vectors, the median is the middle element
        return v[n];
    }
    else {
        // For even-sized vectors, find the average of the two middle elements
        std::vector<float> temp(v.begin(), v.begin() + n + 1);
        std::nth_element(temp.begin(), temp.begin() + n - 1, temp.end());
        return (temp[n - 1] + v[n]) / 2.0f;
    }
}



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

    float sum = std::transform_reduce(std::execution::par, v.begin(), v.end(), 0.0f, std::plus<>(), logTransform);

    return sum / static_cast<float>(positiveCount);
}


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
    _tableModel(this, "Table Model"),
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
    _topNGenesFilter(this, "Top N", 10),
    _geneNamesConnection(this, "Gene Names Connection"),
    _createRowMultiSelectTree(this, "Create Row MultiSelect Tree"),
    _performGeneTableTsneAction(this, "Perform Gene Table TSNE"),
    _tsnePerplexity(this, "TSNE Perplexity"),
    _hiddenShowncolumns(this, "Hidden Shown Columns"),
    _scatterplotReembedColorOption(this, "Reembed Color"),
    _scatterplotEmbeddingPointsUMAPOption(this, "Embedding UMAP Points"),
    _selectedSpeciesVals(this, "Selected Species Vals"),
    _removeRowSelection(this, "Remove Selection"),
    _statusColorAction(this, "Status color"),
    _typeofTopNGenes(this, "N Type"),
    _usePreComputedTSNE(this, "Use Precomputed TSNE")
{
    setSerializationName("CSCGDV:CrossSpeciesComparison Gene Detect Plugin Settings");
    _statusBarActionWidget  = new QStatusBar();
    _tableView = new QTableView();
    _selectionDetailsTable = new QTableView();
    _splitter = new QHBoxLayout();

    _tableView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    _tableView->setSelectionBehavior(QAbstractItemView::SelectRows);
    _tableView->setSelectionMode(QAbstractItemView::SingleSelection);
    _tableView->setEditTriggers(QAbstractItemView::NoEditTriggers);
    _tableView->setAlternatingRowColors(true);
    _tableView->setSortingEnabled(true);
    _tableView->setShowGrid(true);
    _tableView->setGridStyle(Qt::SolidLine);
    _tableView->setHorizontalScrollMode(QAbstractItemView::ScrollPerPixel);
    _tableView->setVerticalScrollMode(QAbstractItemView::ScrollPerPixel);
    _tableView->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _tableView->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _tableView->setCornerButtonEnabled(false);
    _tableView->setWordWrap(false);
    _tableView->setTabKeyNavigation(false);
    _tableView->setAcceptDrops(false);
    _tableView->setDropIndicatorShown(false);
    _tableView->setDragEnabled(false);
    _tableView->setDragDropMode(QAbstractItemView::NoDragDrop);
    _tableView->setDragDropOverwriteMode(false);
    _tableView->setAutoScroll(false);
    _tableView->setAutoScrollMargin(16);
    _tableView->setAutoFillBackground(true);
    _tableView->setFrameShape(QFrame::NoFrame);
    _tableView->setFrameShadow(QFrame::Plain);
    _tableView->setLineWidth(0);
    _tableView->setMidLineWidth(0);
    _tableView->setFocusPolicy(Qt::NoFocus);
    _tableView->setContextMenuPolicy(Qt::NoContextMenu);
    _tableView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    _tableView->setMinimumSize(QSize(0, 0));
    _tableView->setMaximumSize(QSize(16777215, 16777215));
    _tableView->setBaseSize(QSize(0, 0));
    _tableView->setFocusPolicy(Qt::StrongFocus);
    _tableView->setVerticalScrollMode(QAbstractItemView::ScrollPerPixel);

    //only highlight multiple rows if shiuft is pressed
    _tableView->setSelectionBehavior(QAbstractItemView::SelectRows);



    _selectionDetailsTable->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    _selectionDetailsTable->setSelectionBehavior(QAbstractItemView::SelectRows);
    _selectionDetailsTable->setSelectionMode(QAbstractItemView::SingleSelection);
    _selectionDetailsTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    _selectionDetailsTable->setAlternatingRowColors(true);
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

    //only highlight multiple rows if shiuft is pressed
    _selectionDetailsTable->setSelectionBehavior(QAbstractItemView::SelectRows);




    _statusBarActionWidget->setStatusTip("Status");
    _statusBarActionWidget->setFixedHeight(20);
    _statusBarActionWidget->setFixedWidth(100);
    _statusBarActionWidget->setAutoFillBackground(true);
    _statusBarActionWidget->setSizeGripEnabled(false);


    _selectedCellClusterInfoStatusBar = new mv::gui::FlowLayout();


    _tableModel.setSerializationName("CSCGDV:Table Model");
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
    _usePreComputedTSNE.setChecked(false);
    _hiddenShowncolumns.setSerializationName("CSCGDV:Hidden Shown Columns");
    _scatterplotReembedColorOption.setSerializationName("CSCGDV:Scatterplot Reembedding Color Option");
    _scatterplotEmbeddingPointsUMAPOption.setSerializationName("CSCGDV:Scatterplot Embedding UMAP Points Option");
    _typeofTopNGenes.setSerializationName("CSCGDV:Type of Top N Genes");
    _performGeneTableTsneAction.setChecked(false);
    _createRowMultiSelectTree.setDisabled(true);
    _selectedRowIndex.setDisabled(true);
    _selectedRowIndex.setString("");
    _scatterplotReembedColorOption.initialize({"Species","Cluster","Expression"}, "Species");
    _typeofTopNGenes.initialize({"Absolute","Negative","Positive","Mixed"}, "Positive");
   
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
    const auto updateGeneFilteringTrigger = [this]() -> void
        {
            _startComputationTriggerAction.setDisabled(true);
            startCodeTimer("UpdateGeneFilteringTrigger");
            startCodeTimer("Part1");
            

            auto pointsDataset = _mainPointsDataset.getCurrentDataset();
            auto embeddingDataset = _embeddingDataset.getCurrentDataset();
            auto speciesDataset = _speciesNamesDataset.getCurrentDataset();
            auto clusterDataset = _clusterNamesDataset.getCurrentDataset();
            auto referenceTreeDataset = _referenceTreeDataset.getCurrentDataset();
            _selectedSpeciesVals.setString("");
            _geneNamesConnection.setString("");
            bool isValid = false;
            
            QString referenceTreedatasetId = "";
            stopCodeTimer("Part1");
            startCodeTimer("Part2");
            if (!pointsDataset.isValid() || !embeddingDataset.isValid() || !speciesDataset.isValid() || !clusterDataset.isValid() || !referenceTreeDataset.isValid())
            {
                qDebug() << "No datasets selected";
                _startComputationTriggerAction.setDisabled(false);
                return;
            }
            if (pointsDataset->getSelectionIndices().size() <1)
            {
                qDebug() << "No points selected";
                _startComputationTriggerAction.setDisabled(false);
                return;
            }
            if (_selectedPointsTSNEDataset.isValid())
            {
                _selectedPointsTSNEDataset->setSelectionIndices({});
            }
            stopCodeTimer("Part2");
            startCodeTimer("Part3");
            _clusterNameToGeneNameToExpressionValue.clear();
            referenceTreedatasetId = referenceTreeDataset->getId();
            isValid = speciesDataset->getParent() == pointsDataset && clusterDataset->getParent() == pointsDataset && embeddingDataset->getParent() == pointsDataset;
            if (!isValid)
            {
                qDebug() << "Datasets are not valid";
                _startComputationTriggerAction.setDisabled(false);
                return;
            }
            _selectedIndicesFromStorage.clear();
            _selectedIndicesFromStorage = pointsDataset->getSelectionIndices();

            auto embeddingDatasetRaw = mv::data().getDataset<Points>(embeddingDataset->getId());
            auto pointsDatasetRaw = mv::data().getDataset<Points>(pointsDataset->getId());
            auto pointsDatasetallColumnNameList = pointsDatasetRaw->getDimensionNames();
            auto embeddingDatasetallColumnNameList = embeddingDatasetRaw->getDimensionNames();
            stopCodeTimer("Part3");
            startCodeTimer("Part4");
            std::vector<int> embeddingDatasetColumnIndices(embeddingDatasetallColumnNameList.size());
            std::iota(embeddingDatasetColumnIndices.begin(), embeddingDatasetColumnIndices.end(), 0);

            std::vector<int> pointsDatasetallColumnIndices(pointsDatasetallColumnNameList.size());
            std::iota(pointsDatasetallColumnIndices.begin(), pointsDatasetallColumnIndices.end(), 0);
            stopCodeTimer("Part4");
        {
            
                if (_selectedIndicesFromStorage.size() > 0 && embeddingDatasetColumnIndices.size() > 0)
            {
                    startCodeTimer("Part5");
                    auto speciesDatasetRaw = mv::data().getDataset<Clusters>(speciesDataset->getId());
                auto clusterDatasetRaw = mv::data().getDataset<Clusters>(clusterDataset->getId());
                auto clusterDatasetName= clusterDatasetRaw->getGuiName();
                auto clustersValuesAll = clusterDatasetRaw->getClusters();
                auto speciesValuesAll = speciesDatasetRaw->getClusters();

                std::map<QString, std::pair<QColor, std::vector<int>>> selctedClustersMap;
                std::map<QString, std::pair<QColor, std::vector<int>>> selectedSpeciesMap;
                stopCodeTimer("Part5");
                if (!speciesValuesAll.empty() && !clustersValuesAll.empty())
                {
                    startCodeTimer("Part6.1");
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
                    stopCodeTimer("Part6.1");
                    if (_selectedPointsDataset.isValid() && _selectedPointsEmbeddingDataset.isValid() && _tsneDatasetSpeciesColors.isValid() && _tsneDatasetClusterColors.isValid())
                    {
                       startCodeTimer("Part6.2");
                        _tsneDatasetSpeciesColors->getClusters() = QVector<Cluster>();
                        events().notifyDatasetDataChanged(_tsneDatasetSpeciesColors);
                        _tsneDatasetClusterColors->getClusters() = QVector<Cluster>();
                        events().notifyDatasetDataChanged(_tsneDatasetClusterColors);
                        stopCodeTimer("Part6.2");
                        startCodeTimer("Part7");
                        startCodeTimer("Part7.1");
                         int selectedIndicesFromStorageSize = _selectedIndicesFromStorage.size();
                         int pointsDatasetColumnsSize = pointsDatasetallColumnIndices.size();
                         int embeddingDatasetColumnsSize = embeddingDatasetColumnIndices.size();
                         QString datasetIdEmb = _selectedPointsDataset->getId();
                         QString datasetId = _selectedPointsEmbeddingDataset->getId();
                         int dimofDatasetExp = 1;
                         std::vector<QString> dimensionNamesExp = { "Expression" };
                         QString datasetIdExp = _tsneDatasetExpressionColors->getId();
                         stopCodeTimer("Part7.1");
                         startCodeTimer("Part7.2");

                         //first thread start
                        std::vector<float> resultContainerForSelectedPoints(selectedIndicesFromStorageSize* pointsDatasetColumnsSize);
                        pointsDatasetRaw->populateDataForDimensions(resultContainerForSelectedPoints, pointsDatasetallColumnIndices, _selectedIndicesFromStorage);
                        populatePointData(datasetIdEmb, resultContainerForSelectedPoints, selectedIndicesFromStorageSize, pointsDatasetColumnsSize, pointsDatasetallColumnNameList);
                        

                        //second thread start
                        std::vector<float> resultContainerForSelectedEmbeddingPoints(selectedIndicesFromStorageSize* embeddingDatasetColumnsSize);
                        embeddingDatasetRaw->populateDataForDimensions(resultContainerForSelectedEmbeddingPoints, embeddingDatasetColumnIndices, _selectedIndicesFromStorage);
                        populatePointData(datasetId, resultContainerForSelectedEmbeddingPoints, selectedIndicesFromStorageSize, embeddingDatasetColumnsSize, embeddingDatasetallColumnNameList);
                        

                        //third thread start
                        std::vector<float> resultContainerColorPoints(selectedIndicesFromStorageSize, -1.0f);
                        populatePointData(datasetIdExp, resultContainerColorPoints, selectedIndicesFromStorageSize, dimofDatasetExp, dimensionNamesExp);
                        
                        //wait for all threads to finish

                        stopCodeTimer("Part7.2");

                        stopCodeTimer("Part7");
                        startCodeTimer("Part8");
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
                        stopCodeTimer("Part8");
                        startCodeTimer("Part9");
                        mv::plugin::AnalysisPlugin* analysisPlugin; 
                        bool usePreTSNE= _usePreComputedTSNE.isChecked();
                        
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
                        stopCodeTimer("Part9"); 
                        if (!usePreTSNE)
                        {
                           startCodeTimer("Part10");
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
                                    _startComputationTriggerAction.setDisabled(false);
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
                            stopCodeTimer("Part10");
                            }
                        else
                        {
                            startCodeTimer("Part11");
                            auto umapDataset = _scatterplotEmbeddingPointsUMAPOption.getCurrentDataset();

                            if (umapDataset.isValid()) 
                            {
        

                                _selectedPointsTSNEDataset = mv::data().createDerivedDataset<Points>("SelectedPointsTSNEDataset", _selectedPointsEmbeddingDataset, _selectedPointsEmbeddingDataset);

                                mv::events().notifyDatasetAdded(_selectedPointsTSNEDataset);

                                auto umapDatasetRaw = mv::data().getDataset<Points>(umapDataset->getId());
                                auto dimNames = umapDatasetRaw->getDimensionNames();
                                int preComputedEmbeddingColumnsSize = umapDatasetRaw->getNumDimensions();
                                std::vector<float> resultContainerPreComputedUMAP(selectedIndicesFromStorageSize* preComputedEmbeddingColumnsSize);
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


                            stopCodeTimer("Part11");


                           
                        }



                    }
                    else
                    {
                        qDebug() << "Datasets are not valid";
                    }
                    startCodeTimer("Part12");
                    startCodeTimer("Part12.1");
                    QtConcurrent::run([&]() {
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
                                selctedClustersMap[clusterName] = { clusterColor, filteredIndices };
                            }
                        }
                        });

                    stopCodeTimer("Part12.1");
                    startCodeTimer("Part12.2");




                    std::mutex clusterGeneMeanExpressionMapMutex;
                    std::mutex clusterNameToGeneNameToExpressionValueMutex;
                    std::mutex nonselectedClusterMeanExpressionMapMutex;

                    std::vector<std::future<void>> futures;

                    for (auto& species : speciesValuesAll) {
                        futures.push_back(std::async(std::launch::async, [&, species]() { // Fixed Problem 7 by capturing species by value
                            auto speciesIndices = species.getIndices();
                            auto speciesName = species.getName();
                            auto speciesColor = species.getColor();

                            std::vector<int> filteredIndices;
                            for (int i = 0; i < speciesIndices.size(); i++) {
                                int indexVal = findIndex(_selectedIndicesFromStorage, speciesIndices[i]);
                                if (indexVal != -1) {
                                    filteredIndices.push_back(indexVal);
                                }
                            }

                            std::vector<int> commonSelectedIndices;
                            std::vector<int> commonNotSelectedIndices;
                            std::sort(_selectedIndicesFromStorage.begin(), _selectedIndicesFromStorage.end());
                            std::sort(speciesIndices.begin(), speciesIndices.end());
                            std::set_intersection(_selectedIndicesFromStorage.begin(), _selectedIndicesFromStorage.end(), speciesIndices.begin(), speciesIndices.end(), std::back_inserter(commonSelectedIndices));
                            std::set_difference(speciesIndices.begin(), speciesIndices.end(), commonSelectedIndices.begin(), commonSelectedIndices.end(), std::back_inserter(commonNotSelectedIndices));

                            for (int i = 0; i < pointsDatasetallColumnNameList.size(); i++) {
                                auto& geneName = pointsDatasetallColumnNameList[i];
                                std::vector<int> geneIndex = { i }; // Fixed Problem 3 by changing auto to std::vector<int>

                                float fullMean = 0.0f; // Initialize meanValue to prevent using uninitialized memory
                                StatisticsSingle calculateStatisticsShort;
                                StatisticsSingle calculateStatisticsNot;
                                if (!commonSelectedIndices.empty()) {
                                    std::vector<float> resultContainerShort(commonSelectedIndices.size());
                                    pointsDatasetRaw->populateDataForDimensions(resultContainerShort, geneIndex, commonSelectedIndices);

                                    calculateStatisticsShort = calculateStatistics(resultContainerShort);
                                    float selectedMean = calculateStatisticsShort.meanVal;
                                    _selectedSpeciesCellCountMap[speciesName].selectedCellsCount = commonSelectedIndices.size();

                                    // Removed problematic code block related to _clusterGeneMeanExpressionMap
                                }

                                if (!commonNotSelectedIndices.empty()) {
                                    std::vector<float> resultContainerShortNot(commonNotSelectedIndices.size());
                                    pointsDatasetRaw->populateDataForDimensions(resultContainerShortNot, geneIndex, commonNotSelectedIndices);

                                    calculateStatisticsNot = calculateStatistics(resultContainerShortNot);
                                    _selectedSpeciesCellCountMap[speciesName].nonSelectedCellsCount = commonNotSelectedIndices.size();
                                }
                                Statistics combinedValue;
                                // Fixed Problem 1 by removing the problematic if condition
                                combinedValue = combineStatisticsSingle(calculateStatisticsShort, calculateStatisticsNot);

                                {
                                    std::lock_guard<std::mutex> lock(clusterNameToGeneNameToExpressionValueMutex);
                                    _clusterNameToGeneNameToExpressionValue[speciesName][geneName] = combinedValue; // Fixed Problem 2, 4, 5, 6, 8, 9, 10 by correcting the structure and usage of combinedValue
                                }
                            }
                            }));
                    }

                    // Wait for all futures to complete
                    for (auto& future : futures) {
                        future.get();
                    }



                    stopCodeTimer("Part12.2");

                    auto clusterColorDatasetId = _tsneDatasetClusterColors->getId();
                    auto speciesColorDatasetId = _tsneDatasetSpeciesColors->getId();
                    startCodeTimer("Part12.3");
                    populateClusterData(speciesColorDatasetId, selectedSpeciesMap);
                    stopCodeTimer("Part12.3");
                    startCodeTimer("Part12.4"); 
                    populateClusterData(clusterColorDatasetId, selctedClustersMap);
                    stopCodeTimer("Part12.4");
                    stopCodeTimer("Part12");
                    if (_tsneDatasetClusterColors.isValid())
                    {
                        
                        auto clusterValues = _tsneDatasetClusterColors->getClusters();
                        if (!clusterValues.empty())
                        {
                            startCodeTimer("Part13");
                            
                            QLayoutItem* layoutItem;
                            while ((layoutItem = _selectedCellClusterInfoStatusBar->takeAt(0)) != nullptr) {
                                delete layoutItem->widget();
                                delete layoutItem;
                            }
                            
                            // Create a description label
                            auto descriptionLabel = new QLabel("Selected Cell Counts per "+ clusterDatasetName +" :");
                            // Optionally, set a stylesheet for the description label for styling
                            descriptionLabel->setStyleSheet("QLabel { font-weight: bold; padding: 2px; }");
                            // Add the description label to the layout
                            _selectedCellClusterInfoStatusBar->addWidget(descriptionLabel);


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


                        }
                        
                    }

                    startCodeTimer("Part14");
                    QVariant geneListTable = findTopNGenesPerCluster(_clusterNameToGeneNameToExpressionValue, _topNGenesFilter.getValue(), referenceTreedatasetId, 1.0);
                    stopCodeTimer("Part14");
                    if (!geneListTable.isNull())
                    {
                        startCodeTimer("Part15");
                        //_filteredGeneNamesVariant.setVariant(geneListTable);
                        _tableModel.setVariant(geneListTable);
                        stopCodeTimer("Part15");

                    }
                    else
                    {
                        qDebug() << "QVariant empty";
                    }



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
        }
        stopCodeTimer("UpdateGeneFilteringTrigger");
        _startComputationTriggerAction.setDisabled(false);
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

 
        computeGeneMeanExpressionMap();
        if (_mainPointsDataset.getCurrentDataset().isValid())
        {
            
            auto fullDataset=mv::data().getDataset<Points>(_mainPointsDataset.getCurrentDataset().getDatasetId());
            auto dimensions = fullDataset->getNumDimensions();
            if (dimensions>0)
            {
                _topNGenesFilter.setMinimum(0);
                _topNGenesFilter.setMaximum(dimensions);

                _topNGenesFilter.setValue(std::min(10, static_cast<int>(dimensions)));

            }
            else

            {
                _topNGenesFilter.setMinimum(0);
                _topNGenesFilter.setMaximum(0);
                _topNGenesFilter.setValue(0);
            }

        }
        else
        {
            _topNGenesFilter.setMinimum(0);
            _topNGenesFilter.setMaximum(0);
            _topNGenesFilter.setValue(0);
            
        }
        
 };

    connect(&_mainPointsDataset, &DatasetPickerAction::currentIndexChanged, this, updateMainPointsDataset);

    const auto updateSpeciesNameDataset = [this]() -> void {

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

                    }

                }

            }

        }
        
        
        computeGeneMeanExpressionMap();
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
        // Assuming the context is to modify the QStatusBar _statusBarActionWidget based on the string value
        auto string = _statusColorAction.getString();
        QString labelText = "";
        QString backgroundColor = "none";
        if (string == "C") {
            labelText = "Updated";
            backgroundColor = "#28a745"; // Green
        }
        else if (string == "M") {
            labelText = "Pending";
            backgroundColor = "#ffc107"; // Gold
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
    const auto updateTopGenesSlider = [this]() -> void {
        _statusColorAction.setString("M");

        };
    connect(&_topNGenesFilter, &IntegralAction::valueChanged, this, updateTopGenesSlider);

    _statusColorAction.setString("M");
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
            QFuture<void> future = QtConcurrent::map(speciesclusters.begin(), speciesclusters.end(), [&](const auto& species) {
                auto speciesIndices = species.getIndices();
                auto speciesName = species.getName();
                for (int i = 0; i < mainPointDimensionNames.size(); i++) {
                    auto& geneName = mainPointDimensionNames[i];
                    auto geneIndex = { i };
                    std::vector<float> resultContainerFull(speciesIndices.size());
                    mainPointDatasetFull->populateDataForDimensions(resultContainerFull, geneIndex, speciesIndices);
                    float fullMean = calculateMean(resultContainerFull);
                    _clusterGeneMeanExpressionMap[speciesName][geneName] = fullMean;
                }
                });
            future.waitForFinished();
        }
    }

}

QVariant SettingsAction::findTopNGenesPerCluster(const std::map<QString, std::map<QString, Statistics>>& map, int n, QString datasetId, float treeSimilarityScore) {

    if (map.empty() || n <= 0) {
        return QVariant();
    }
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

    for (const auto& outerPair : map) {
        
        auto speciesName=outerPair.first;        
        std::vector<std::pair<QString, float>> geneExpressionVec;
        geneExpressionVec.reserve(outerPair.second.size());
        for (const auto& innerPair : outerPair.second) {

            auto geneName = innerPair.first;
            auto differenceMeanValue = innerPair.second.meanSelected-innerPair.second.meanNonSelected;
            geneExpressionVec.push_back(std::make_pair(geneName, differenceMeanValue));

        }
        // Sort the geneExpressionVec based on the mean value from highest to lowest
        std::sort(geneExpressionVec.begin(), geneExpressionVec.end(), [](const auto& a, const auto& b) {
            return a.second > b.second;
            });

        // Select the top n genes based on sorted list. If Positive then top n genes with highest mean value, if Negative then top n genes with lowest mean value, if Mixed then top n/2 genes with highest mean value and n/2 genes with lowest mean value, if Absolute then top n genes with highest absolute mean value.
        switch (option) {
        case SelectionOption::AbsoluteTopN: {
            std::sort(geneExpressionVec.begin(), geneExpressionVec.end(), [](const auto& a, const auto& b) {
                return std::abs(a.second) > std::abs(b.second);
                });
            for (int i = 0; i < std::min(n, static_cast<int>(geneExpressionVec.size())); ++i) {
                geneList.insert(geneExpressionVec[i].first);
                geneAppearanceCounter[geneExpressionVec[i].first].push_back(outerPair.first);
            }
            break;
        }
        case SelectionOption::PositiveTopN: {
            for (int i = 0; i < std::min(n, static_cast<int>(geneExpressionVec.size())); ++i) {
                geneList.insert(geneExpressionVec[i].first);
                geneAppearanceCounter[geneExpressionVec[i].first].push_back(outerPair.first);
            }
            break;
        }
        case SelectionOption::NegativeTopN: {
            std::reverse(geneExpressionVec.begin(), geneExpressionVec.end()); // Reverse to get the lowest values
            for (int i = 0; i < std::min(n, static_cast<int>(geneExpressionVec.size())); ++i) {
                geneList.insert(geneExpressionVec[i].first);
                geneAppearanceCounter[geneExpressionVec[i].first].push_back(outerPair.first);
            }
            break;
        }
        case SelectionOption::MixedTopN: {
            int halfN = n / 2;
            for (int i = 0; i < halfN; ++i) { // Top halfN genes with highest mean value
                geneList.insert(geneExpressionVec[i].first);
                geneAppearanceCounter[geneExpressionVec[i].first].push_back(outerPair.first);
            }
            if (n % 2 != 0) { // If n is odd, add one more from the top half
                geneList.insert(geneExpressionVec[halfN].first);
                geneAppearanceCounter[geneExpressionVec[halfN].first].push_back(outerPair.first);
                halfN++;
            }
            for (int i = geneExpressionVec.size() - halfN; i < geneExpressionVec.size(); ++i) { // Bottom halfN genes with lowest mean value
                geneList.insert(geneExpressionVec[i].first);
                geneAppearanceCounter[geneExpressionVec[i].first].push_back(outerPair.first);
            }
            break;
        }
        }
    }

    return createModelFromData(geneList, map, datasetId, treeSimilarityScore, geneAppearanceCounter, n);
}

QVariant SettingsAction::createModelFromData(const QSet<QString>& returnGeneList, const std::map<QString, std::map<QString, Statistics>>& map, const QString& treeDatasetId, const float& treeSimilarityScore, const std::map<QString, std::vector<QString>>& geneCounter, const int& n) {

    if (returnGeneList.isEmpty() || map.empty()) {
        return QVariant();
    }

    const std::unordered_map<std::string, int> clusteringTypeMap = {
    {"Complete", HCLUST_METHOD_COMPLETE},
    {"Average", HCLUST_METHOD_AVERAGE},
    {"Median", HCLUST_METHOD_MEDIAN},
    {"Single", HCLUST_METHOD_SINGLE} // Added "Single" to the map for consistency
    };

    QStandardItemModel* model = new QStandardItemModel();
    int numOfSpecies = map.size();
    _initColumnNames.clear();
    _initColumnNames = { "ID", /*"Newick tree", "Similarity with Reference Tree",*/ "Mean \nExpression", "Species \nAppearance", "Gene Appearance Species Names" ,"Statistics"};
    model->setHorizontalHeaderLabels(_initColumnNames);

    QStringList headers = _initColumnNames;
    _hiddenShowncolumns.setSelectedOptions({});
    _hiddenShowncolumns.setOptions({});
    _hiddenShowncolumns.setOptions(headers);

    QStringList selectedHeaders = { headers[0], headers[1], headers[2] };//{ headers[0], headers[2], headers[3], headers[4] };
    _hiddenShowncolumns.setSelectedOptions(selectedHeaders);


    std::map<QString, std::pair<QString, std::map<QString, Statistics>>> newickTrees;

    std::string clusteringTypecurrentText = "Single";  // "Single", "Complete", "Average", "Median"
    int opt_method = clusteringTypeMap.at(clusteringTypecurrentText); 

    for (const auto& gene : returnGeneList) {
        QList<QStandardItem*> row;
        std::vector<float> numbers;
        std::map<QString, Statistics> statisticsValuesForSpeciesMap;

        for (const auto& outerPair : map) {
            QString speciesName = outerPair.first;

            const std::map<QString, Statistics>& innerMap = outerPair.second;
            auto it = innerMap.find(gene);
            float value = 0.0;
           // QString typeCal = "Mean";
            //if (typeCal == "Median")
            //{
                //value = it->second.median;
            //}
            //else //Mean
            //{
                value = it->second.meanSelected;
            //}
            numbers.push_back(value);
            statisticsValuesForSpeciesMap[speciesName] = it->second;
        }

        auto numOfLeaves = static_cast<int>(numOfSpecies);
        std::unique_ptr<double[]> distmat(condensedDistanceMatrix(numbers));
        std::unique_ptr<int[]> merge(new int[2 * (numOfLeaves - 1)]);
        std::unique_ptr<double[]> height(new double[numOfLeaves - 1]);
        hclust_fast(numOfLeaves, distmat.get(), opt_method, merge.get(), height.get());
        std::string newick = mergeToNewick(merge.get(), numOfLeaves);
        newickTrees.insert({ gene, {QString::fromStdString(newick), statisticsValuesForSpeciesMap} });

        row.push_back(new QStandardItem(gene)); //0 ID
        //row.push_back(new QStandardItem(""));  //1 Newick tree
        //row.push_back(new QStandardItem(QString::number(-1.0)));  //2 Similarity with Reference Tree

        float meanV = 0.0;
        //iterate statisticsValuesForSpeciesMap
        for (const auto& pair : statisticsValuesForSpeciesMap) {
            meanV += pair.second.meanSelected;
        }
        meanV= meanV / statisticsValuesForSpeciesMap.size();
        row.push_back(new QStandardItem(QString::number(meanV)));//3 Mean Expression

        auto it = geneCounter.find(gene);
        QString speciesGeneAppearancesComb;
        int count = 0;
        if (it != geneCounter.end()) {
            count = it->second.size();
            speciesGeneAppearancesComb = QStringList(it->second.begin(), it->second.end()).join(";");
        }
        else
        {
            qDebug() << "Gene not found in geneCounter";
        }
        row.push_back(new QStandardItem(QString::number(count))); // 4 Gene Appearances /" + QString::number(numOfSpecies) + " Species"
        row.push_back(new QStandardItem(speciesGeneAppearancesComb)); // 5 Gene Appearance Species Names

        QString formattedStatistics;
        for (const auto& pair : statisticsValuesForSpeciesMap) {
            const auto& species = pair.first;
            const auto& stats = pair.second;
            formattedStatistics += QString("Species: %1, MeanSelected: %2, MedianSelected: %3, ModeSelected: %4, RangeSelected: %5, CountSelected: %6, MeanNotSelected: %7, MedianNotSelected: %8, ModeNotSelected: %9, RangeNotSelected: %10, CountNotSelected: %11;\n")
                .arg(species)
                .arg(stats.meanSelected, 0, 'f', 2)
                .arg(stats.medianSelected, 0, 'f', 2)
                .arg(stats.modeSelected, 0, 'f', 2)
                .arg(stats.rangeSelected, 0, 'f', 2)
                .arg(stats.countSelected)
                .arg(stats.meanNonSelected, 0, 'f', 2)
                .arg(stats.medianNonSelected, 0, 'f', 2)
                .arg(stats.modeNonSelected, 0, 'f', 2)
                .arg(stats.rangeNonSelected, 0, 'f', 2)
                .arg(stats.countNonSelected);
        }

        row.push_back(new QStandardItem(formattedStatistics)); //6 Statistics

        model->appendRow(row);
    }


    /*
    std::string targetNewick = "";
    QStringList fullTreeNames;

    std::vector<QString> leafnames;
    leafnames.reserve(map.size());
    for (const auto& outerPair : map) {
        leafnames.push_back(outerPair.first);
    }

    std::map<QString, std::pair<QString, float>> treeSimilarities;
    if (!treeDatasetId.isEmpty()) 
    {
        auto fullTreeData = mv::data().getDataset<CrossSpeciesComparisonTree>(treeDatasetId);
        if (fullTreeData.isValid())
        {
            fullTreeNames = fullTreeData->getTreeLeafNames();
            auto treeData = fullTreeData->getTreeData();
            auto jsonTree = nlohmann::json::parse(QJsonDocument(treeData).toJson(QJsonDocument::Compact).toStdString());
            targetNewick = jsonToNewick(jsonTree, leafnames);
            targetNewick += ";";  // End of Newick string
        }
    }



    if (fullTreeNames.size() > 0 && leafnames.size() > 0 && targetNewick != "")
    {

        QStringList copyleafNames = QStringList(leafnames.begin(), leafnames.end());

        if (areSameIgnoreOrder(fullTreeNames, copyleafNames)) {
  
            for (auto& pair : newickTrees) {//need to change

                std::string modifiedNewick = pair.second.first.toStdString();
                std::map <QString, Statistics> speciesStatisticsMaps = pair.second.second;

                const char* string1 = targetNewick.c_str();
                const char* string2 = modifiedNewick.c_str();
                Tree t1;
                Tree t2;

                // Fix for Problem 1: Check the return value of freopen
                if (freopen("CON", "r", stdin) == nullptr) {
                    std::cerr << "Failed to reopen stdin from CON." << std::endl;
                }

                FILE* file1 = nullptr;
                FILE* file2 = nullptr;

                // Fix for Problem 2: Ensure file1 is not nullptr before using fputs
                if (fopen_s(&file1, "file1.txt", "w") == 0 && file1 != nullptr) {
                    fputs(string1, file1); // Ensure string1 is a std::string or converted to C-string
                    fclose(file1);
                }
                else {
                    std::cerr << "Failed to open file1.txt for writing." << std::endl;
                }

                // Fix for Problem 3: Ensure file2 is not nullptr before using fputs
                if (fopen_s(&file2, "file2.txt", "w") == 0 && file2 != nullptr) {
                    fputs(string2, file2); // Ensure string2 is a std::string or converted to C-string
                    fclose(file2);
                }
                else {
                    std::cerr << "Failed to open file2.txt for writing." << std::endl;
                }

                // Fix for Problem 4: Check the return value of freopen
                if (freopen("file1.txt", "r", stdin) == nullptr) {
                    std::cerr << "Failed to reopen stdin from file1.txt." << std::endl;
                }
                else {
                    t1.CreateTree();
                }

                // Fix for Problem 5: Check the return value of freopen
                if (freopen("file2.txt", "r", stdin) == nullptr) {
                    std::cerr << "Failed to reopen stdin from file2.txt." << std::endl;
                }
                else {
                    t2.CreateTree();
                }

                int sim = Calculate(&t1, &t2);


                float similarity = 1.0 - static_cast<float>(sim) / static_cast<float>(numOfSpecies); 
                std::pair<QString, float>  temp;
                temp.first = createJsonTreeFromNewick(QString::fromStdString(modifiedNewick), leafnames, speciesStatisticsMaps);
                temp.second = similarity;
                treeSimilarities.insert(std::make_pair(pair.first, temp));

            }

        }
    }

    for (int i = 0; i < model->rowCount(); ++i) {
        QString gene = model->item(i, 0)->text();

        auto it = treeSimilarities.find(gene);
        if (it != treeSimilarities.end()) {
            const auto& [newick, similarity] = it->second;

            model->item(i, 1)->setText(newick);

            auto item2 = model->item(i, 2);
            item2->setData(similarity, Qt::DisplayRole);
            item2->setData(similarity, Qt::UserRole);
        }
    }
    */


    return QVariant::fromValue(model);

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
    _tableModel.fromParentVariantMap(variantMap);
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
    _scatterplotReembedColorOption.fromParentVariantMap(variantMap);
    _scatterplotEmbeddingPointsUMAPOption.fromParentVariantMap(variantMap);
    _selectedSpeciesVals.fromParentVariantMap(variantMap);
    _removeRowSelection.fromParentVariantMap(variantMap);
    _statusColorAction.fromParentVariantMap(variantMap);
    _typeofTopNGenes.fromParentVariantMap(variantMap);
    _usePreComputedTSNE.fromParentVariantMap(variantMap);
 

}

QVariantMap SettingsAction::toVariantMap() const
{
    QVariantMap variantMap = WidgetAction::toVariantMap();

    _geneNamesConnection.insertIntoVariantMap(variantMap);
    _createRowMultiSelectTree.insertIntoVariantMap(variantMap);
    _tableModel.insertIntoVariantMap(variantMap);
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
    _scatterplotReembedColorOption.insertIntoVariantMap(variantMap);
    _scatterplotEmbeddingPointsUMAPOption.insertIntoVariantMap(variantMap);
    _selectedSpeciesVals.insertIntoVariantMap(variantMap);
    _removeRowSelection.insertIntoVariantMap(variantMap);
    _statusColorAction.insertIntoVariantMap(variantMap);
    _typeofTopNGenes.insertIntoVariantMap(variantMap);
    _usePreComputedTSNE.insertIntoVariantMap(variantMap);
    return variantMap;
}
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
#include <vector>
using namespace mv;
using namespace mv::gui;

float calculateMean(const std::vector<float>& v) {
    if (v.empty())
        return 0.0f;

    float sum = std::reduce(std::execution::par, v.begin(), v.end(), 0.0f);
    return sum / static_cast<float>(v.size());
}

//struct Statistics {
//    float mean;
//    float variance;
//    float stdDeviation;
//};
//Statistics calculateStatistics(const std::vector<float>& numbers) {
//    if (numbers.empty()) {
//        return { 0.0, 0.0, 0.0 };
//    }
//
//    float sum = std::accumulate(numbers.begin(), numbers.end(), 0.0);
//    float mean = sum / numbers.size();
//    float sq_sum = std::inner_product(numbers.begin(), numbers.end(), numbers.begin(), 0.0);
//    float variance = sq_sum / numbers.size() - mean * mean;
//    float stdDeviation = std::sqrt(variance);
//
//    return { mean, variance, stdDeviation };
//}
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
    _speciesNamesDataset(this, "Species Names Dataset"),
    _clusterNamesDataset(this, "Cluster Names Dataset"),
    //_calculationReferenceCluster(this, "Calculation Reference Cluster"),
    _filteredGeneNamesVariant(this, "Filtered Gene Names"),
    _topNGenesFilter(this, "Top N", 10),
    _geneNamesConnection(this, "Gene Names Connection"),
    _createRowMultiSelectTree(this, "Create Row MultiSelect Tree"),
    _performGeneTableTsneAction(this, "Perform Gene Table TSNE"),
    _tsnePerplexity(this, "TSNE Perplexity"),
    _hiddenShowncolumns(this, "Hidden Shown Columns"),
    _scatterplotReembedColorOption(this, "Reembed Color"),
    _scatterplotEmbeddingColorOption(this, "Embed Color"),
    _scatterplotEmbeddingPointsUMAPOption(this, "Embedding UMAP Points"),
    _selectedSpeciesVals(this, "Selected Species Vals"),
    _removeRowSelection(this, "Remove Selection"),
    _statusColorAction(this, "Status color"),
    _typeofTopNGenes(this, "N Type"),
    _usePreComputedTSNE(this, "Use Precomputed TSNE"),
    _selectedCellClusterInfoBox(this, "Selected Cell Cluster Info")
{
    setSerializationName("CSCGDV:CrossSpeciesComparison Gene Detect Plugin Settings");
    _statusBarActionWidget  = new QStatusBar();
    _statusBarActionWidget->setStatusTip("Status");
    _statusBarActionWidget->setFixedHeight(20);
    _statusBarActionWidget->setFixedWidth(100);
    _statusBarActionWidget->setAutoFillBackground(true);
    _statusBarActionWidget->setSizeGripEnabled(false);


    _selectedCellClusterInfoStatusBar = new QStatusBar();
    _selectedCellClusterInfoStatusBar->setStatusTip("Status");
    _selectedCellClusterInfoStatusBar->setMinimumHeight(40);
    _selectedCellClusterInfoStatusBar->setMinimumWidth(100);
    _selectedCellClusterInfoStatusBar->setAutoFillBackground(true);
    _selectedCellClusterInfoStatusBar->setSizeGripEnabled(false);


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
    _selectedCellClusterInfoBox.setSerializationName("CSCGDV:Selected Cell Cluster Info");
    _selectedGene.setDisabled(true);
    _selectedGene.setString("");
    _startComputationTriggerAction.setSerializationName("CSCGDV:Start Computation");
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
    _scatterplotEmbeddingColorOption.setSerializationName("CSCGDV:Scatterplot Embedding Color Option"); 
    _scatterplotEmbeddingPointsUMAPOption.setSerializationName("CSCGDV:Scatterplot Embedding UMAP Points Option");
    _typeofTopNGenes.setSerializationName("CSCGDV:Type of Top N Genes");
    _performGeneTableTsneAction.setChecked(false);
    _createRowMultiSelectTree.setDisabled(true);
    _selectedRowIndex.setDisabled(true);
    _selectedRowIndex.setString("");
    _scatterplotReembedColorOption.initialize({"Species","Cluster","Expression"}, "Species");
    _typeofTopNGenes.initialize({"Absolute","Negative","Positive","Mixed"}, "Absolute");
   
    _scatterplotEmbeddingColorOption.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == ClusterType;
        });
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
            
            auto pointsDataset = _mainPointsDataset.getCurrentDataset();
            auto embeddingDataset = _embeddingDataset.getCurrentDataset();
            auto speciesDataset = _speciesNamesDataset.getCurrentDataset();
            auto clusterDataset = _clusterNamesDataset.getCurrentDataset();
            auto referenceTreeDataset = _referenceTreeDataset.getCurrentDataset();
            _selectedSpeciesVals.setString("");
            _geneNamesConnection.setString("");
            _selectedCellClusterInfoBox.setString("");
            bool isValid = false;
            QString referenceTreedatasetId = "";

            if (!pointsDataset.isValid() || !embeddingDataset.isValid() || !speciesDataset.isValid() || !clusterDataset.isValid() || !referenceTreeDataset.isValid())
            {
                qDebug() << "No datasets selected";
                return;
            }
            if (pointsDataset->getSelectionIndices().size() <1)
            {
                qDebug() << "No points selected";
                return;
            }

            

            if (_selectedPointsTSNEDataset.isValid())
            {
                _selectedPointsTSNEDataset->setSelectionIndices({});
            }
            _clusterNameToGeneNameToExpressionValue.clear();
            referenceTreedatasetId = referenceTreeDataset->getId();
            isValid = speciesDataset->getParent() == pointsDataset && clusterDataset->getParent() == pointsDataset && embeddingDataset->getParent() == pointsDataset;
            if (!isValid)
            {
                qDebug() << "Datasets are not valid";
                return;
            }
            _selectedIndicesFromStorage.clear();
            _selectedIndicesFromStorage = pointsDataset->getSelectionIndices();

            auto embeddingDatasetRaw = mv::data().getDataset<Points>(embeddingDataset->getId());
            auto pointsDatasetRaw = mv::data().getDataset<Points>(pointsDataset->getId());
            auto pointsDatasetallColumnNameList = pointsDatasetRaw->getDimensionNames();
            auto embeddingDatasetallColumnNameList = embeddingDatasetRaw->getDimensionNames();

            std::vector<int> embeddingDatasetColumnIndices(embeddingDatasetallColumnNameList.size());
            std::iota(embeddingDatasetColumnIndices.begin(), embeddingDatasetColumnIndices.end(), 0);

            std::vector<int> pointsDatasetallColumnIndices(pointsDatasetallColumnNameList.size());
            std::iota(pointsDatasetallColumnIndices.begin(), pointsDatasetallColumnIndices.end(), 0);

        {
            if (_selectedIndicesFromStorage.size() > 0 && embeddingDatasetColumnIndices.size() > 0)
            {
                auto speciesDatasetRaw = mv::data().getDataset<Clusters>(speciesDataset->getId());
                auto clusterDatasetRaw = mv::data().getDataset<Clusters>(clusterDataset->getId());

                auto clustersValuesAll = clusterDatasetRaw->getClusters();
                auto speciesValuesAll = speciesDatasetRaw->getClusters();

                std::map<QString, std::pair<QColor, std::vector<int>>> selctedClustersMap;
                std::map<QString, std::pair<QColor, std::vector<int>>> selectedSpeciesMap;

                if (!speciesValuesAll.empty() && !clustersValuesAll.empty())
                {

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

                    if (_selectedPointsDataset.isValid() && _selectedPointsEmbeddingDataset.isValid() && _tsneDatasetSpeciesColors.isValid() && _tsneDatasetClusterColors.isValid())
                    {
                        _tsneDatasetSpeciesColors->getClusters() = QVector<Cluster>();
                        events().notifyDatasetDataChanged(_tsneDatasetSpeciesColors);
                        _tsneDatasetClusterColors->getClusters() = QVector<Cluster>();
                        events().notifyDatasetDataChanged(_tsneDatasetClusterColors);

                         int selectedIndicesFromStorageSize = _selectedIndicesFromStorage.size();
                         int pointsDatasetColumnsSize = pointsDatasetallColumnIndices.size();
                         int embeddingDatasetColumnsSize = embeddingDatasetColumnIndices.size();

                        // Reserve space for vectors to avoid reallocations
                        std::vector<float> resultContainerForSelectedPoints(selectedIndicesFromStorageSize* pointsDatasetColumnsSize);
                        pointsDatasetRaw->populateDataForDimensions(resultContainerForSelectedPoints, pointsDatasetallColumnIndices, _selectedIndicesFromStorage);

                        QString datasetIdEmb = _selectedPointsDataset->getId();
                        populatePointData(datasetIdEmb, resultContainerForSelectedPoints, selectedIndicesFromStorageSize, pointsDatasetColumnsSize, pointsDatasetallColumnNameList);

                        std::vector<float> resultContainerForSelectedEmbeddingPoints(selectedIndicesFromStorageSize* embeddingDatasetColumnsSize);
                        embeddingDatasetRaw->populateDataForDimensions(resultContainerForSelectedEmbeddingPoints, embeddingDatasetColumnIndices, _selectedIndicesFromStorage);

                        QString datasetId = _selectedPointsEmbeddingDataset->getId();
                        populatePointData(datasetId, resultContainerForSelectedEmbeddingPoints, selectedIndicesFromStorageSize, embeddingDatasetColumnsSize, embeddingDatasetallColumnNameList);

                        std::vector<float> resultContainerColorPoints(selectedIndicesFromStorageSize, -1.0f); 

                        QString datasetIdExp = _tsneDatasetExpressionColors->getId();
                        int dimofDatasetExp = 1;
                        std::vector<QString> dimensionNamesExp = { "Expression" };

                        populatePointData(datasetIdExp, resultContainerColorPoints, selectedIndicesFromStorageSize, dimofDatasetExp, dimensionNamesExp);

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
                        
                        if (!usePreTSNE)
                        {
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
                            }
                        else
                        {
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





                           
                        }



                    }
                    else
                    {
                        qDebug() << "Datasets are not valid";
                    }
                    for (auto& clusters : clustersValuesAll)
                    {
                        auto clusterIndices = clusters.getIndices();
                        auto clusterName = clusters.getName();
                        auto clusterColor = clusters.getColor();
                        std::vector<int> filteredIndices;
                        for (int i = 0; i < clusterIndices.size(); i++)
                        {

                            int indexVal = findIndex(_selectedIndicesFromStorage, clusterIndices[i]);
                            if (indexVal != -1)
                            {
                                filteredIndices.push_back(indexVal);
                            }

                        }
                        selctedClustersMap[clusterName] = { clusterColor, filteredIndices };
                    }


                    for (auto& species : speciesValuesAll)
                    {
                        auto speciesIndices = species.getIndices();
                        auto speciesName = species.getName();
                        auto speciesColor = species.getColor();


                        std::vector<int> filteredIndices;
                        for (int i = 0; i < speciesIndices.size(); i++)
                            {
                                int indexVal = findIndex(_selectedIndicesFromStorage, speciesIndices[i]);
                                if (indexVal != -1)
                                {
                                    filteredIndices.push_back(indexVal);
                                }   
                            }

                        selectedSpeciesMap[speciesName] = { speciesColor, filteredIndices };



                        std::vector<int> commonSelectedIndices;

                        std::sort(_selectedIndicesFromStorage.begin(), _selectedIndicesFromStorage.end());
                        std::sort(speciesIndices.begin(), speciesIndices.end());
                        std::set_intersection(_selectedIndicesFromStorage.begin(), _selectedIndicesFromStorage.end(), speciesIndices.begin(), speciesIndices.end(), std::back_inserter(commonSelectedIndices));






                        for (int i = 0; i < pointsDatasetallColumnNameList.size(); i++) {
                            auto& geneName = pointsDatasetallColumnNameList[i];
                            auto geneIndex = { i };
                            
                            float fullMean;
                            float meanValue;
                            if (!commonSelectedIndices.empty()) {
                                std::vector<float> resultContainerShort(commonSelectedIndices.size());
                                pointsDatasetRaw->populateDataForDimensions(resultContainerShort, geneIndex, commonSelectedIndices);
                                float shortMean = calculateMean(resultContainerShort);
                                
                                
                                auto speciesIter = _clusterGeneMeanExpressionMap.find(speciesName);
                                if (speciesIter == _clusterGeneMeanExpressionMap.end() || speciesIter->second.find(geneName) == speciesIter->second.end()) {

                                    std::vector<float> resultContainerFull(speciesIndices.size());
                                    pointsDatasetRaw->populateDataForDimensions(resultContainerFull, geneIndex, speciesIndices);
                                    fullMean = calculateMean(resultContainerFull);
                                    // Insert the calculated mean into the map
                                    _clusterGeneMeanExpressionMap[speciesName][geneName] = fullMean;
                                    
                                }
                                else
                                {
                                    fullMean= _clusterGeneMeanExpressionMap[speciesName][geneName];
                                }
                                meanValue = fullMean - shortMean;
                            }

                            _clusterNameToGeneNameToExpressionValue[speciesName][geneName] = meanValue;
                        }



                    }


                    auto clusterColorDatasetId = _tsneDatasetClusterColors->getId();
                    auto speciesColorDatasetId = _tsneDatasetSpeciesColors->getId();

                    populateClusterData(speciesColorDatasetId, selectedSpeciesMap);
                    populateClusterData(clusterColorDatasetId, selctedClustersMap);

                    if (_tsneDatasetClusterColors.isValid())
                    {
                        auto clusterValues = _tsneDatasetClusterColors->getClusters();
                        if (!clusterValues.empty())
                        {
                            QString selectedClusterInfo = "<html><head/><body><p>Cell counts: <br>";
                            int clusterCounter = 0;
                            for (auto cluster : clusterValues) {
                                auto clusterName = cluster.getName();
                                auto clusterIndicesSize = cluster.getIndices().size();
                                auto clusterColor = cluster.getColor().name(); // Assuming getColor() returns a QColor

                                // Use the background-color style for the span to set the background of the text
                                selectedClusterInfo += QString(" <span style=\"background-color:%3;\">%1: %2;</span> ").arg(clusterName).arg(clusterIndicesSize).arg(clusterColor);
                                ++clusterCounter;
                                if (clusterCounter % 7 == 0) {
                                    selectedClusterInfo += "<br>"; // Add a new line after every 5 clusters
                                }
                            }
                            if (clusterCounter % 5 != 0) { // If the last line doesn't end with a "; "
                                selectedClusterInfo.chop(2); // Remove the last "; " for formatting
                            }
                            selectedClusterInfo += "</p></body></html>";
                            _selectedCellClusterInfoBox.setString(selectedClusterInfo);
                        }

                    }


                    QVariant geneListTable = findTopNGenesPerCluster(_clusterNameToGeneNameToExpressionValue, _topNGenesFilter.getValue(), referenceTreedatasetId, 1.0);

                    if (!geneListTable.isNull())
                    {
                        //_filteredGeneNamesVariant.setVariant(geneListTable);
                        _tableModel.setVariant(geneListTable);

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

    const auto updateSelectedCellClusterInfoBox = [this]() -> void {

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




        };
    connect(&_selectedCellClusterInfoBox, &StringAction::stringChanged, this, updateSelectedCellClusterInfoBox);
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

QVariant SettingsAction::findTopNGenesPerCluster(const std::map<QString, std::map<QString, float>>& map, int n, QString datasetId, float treeSimilarityScore) {

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
    QStringList returnGeneList;
    std::map<QString, std::vector<QString>> geneAppearanceCounter;

    for (const auto& outerPair : map) {
        std::vector<std::pair<QString, float>> geneExpressionVec;
        // Reserve space to avoid reallocations
        geneExpressionVec.reserve(outerPair.second.size());

        // Filter based on the selection option to avoid unnecessary sorting
        if (option == SelectionOption::NegativeTopN) {
            std::copy_if(outerPair.second.begin(), outerPair.second.end(), std::back_inserter(geneExpressionVec),
                [](const auto& pair) { return pair.second < 0; });
        }
        else {
            geneExpressionVec.assign(outerPair.second.begin(), outerPair.second.end());
        }

        // Lambda for sorting based on the selection option
        auto sortLambda = [&option](const auto& a, const auto& b) {
            if (option == SelectionOption::AbsoluteTopN) {
                return std::abs(a.second) > std::abs(b.second);
            }
            else if (option == SelectionOption::PositiveTopN || option == SelectionOption::MixedTopN) {
                return a.second > b.second;
            }
            else { // NegativeTopN
                return a.second < b.second;
            }
            };

        // Sort once for all options
        std::sort(geneExpressionVec.begin(), geneExpressionVec.end(), sortLambda);

        // Determine the limit for selection based on the option
        int limit = n;
        if (option == SelectionOption::MixedTopN) {
            limit = std::min(n / 2, static_cast<int>(geneExpressionVec.size()));
        }
        else {
            limit = std::min(n, static_cast<int>(geneExpressionVec.size()));
        }

        // Select top elements based on the selection option
        for (int i = 0; i < limit; ++i) {
            const auto& gene = geneExpressionVec[i].first;
            geneList.insert(gene);
            geneAppearanceCounter[gene].push_back(outerPair.first);
        }

        // Additional loop for MixedTopN to select negative values
        if (option == SelectionOption::MixedTopN) {
            int negativeStart = std::max(static_cast<int>(geneExpressionVec.size()) - n / 2, 0);
            for (int i = negativeStart; i < geneExpressionVec.size(); ++i) {
                const auto& gene = geneExpressionVec[i].first;
                geneList.insert(gene);
                geneAppearanceCounter[gene].push_back(outerPair.first);
            }
        }
    }


    returnGeneList = QStringList(geneList.begin(), geneList.end());




    QVariant returnValue = createModelFromData(returnGeneList, map, datasetId, treeSimilarityScore, geneAppearanceCounter, n);

    return returnValue;
}

QVariant SettingsAction::createModelFromData(const QStringList& returnGeneList, const std::map<QString, std::map<QString, float>>& map, const QString& treeDatasetId, const float& treeSimilarityScore, const std::map<QString, std::vector<QString>>& geneCounter, const int& n) {

    if (returnGeneList.isEmpty() || map.empty()) {
        return QVariant();
    }

    QStandardItemModel* model = new QStandardItemModel();
    int numOfSpecies = map.size();
    _initColumnNames.clear();
    _initColumnNames = { "ID", "Newick tree", "Similarity with Reference Tree", "Mean Differential Expression", "Gene Appearances /" + QString::number(numOfSpecies) + " Species", "Gene Appearance Species Names" };
    model->setHorizontalHeaderLabels(_initColumnNames);

    for (auto it = map.cbegin(); it != map.cend(); ++it) {
        QString headerTitle = it->first;
        model->setHorizontalHeaderItem(_initColumnNames.size() + std::distance(map.cbegin(), it), new QStandardItem(headerTitle));
    }

    QStringList headers = _initColumnNames;
    headers.reserve(_initColumnNames.size() + map.size());
    for (auto it = map.cbegin(); it != map.cend(); ++it) {
        headers.push_back(it->first);
    }
    _hiddenShowncolumns.setOptions(headers);

    QStringList selectedHeaders = { headers[0], headers[2], headers[3], headers[4] };
    _hiddenShowncolumns.setSelectedOptions(selectedHeaders);


    std::map<QString, std::pair<QString, std::map<QString, float>>> newickTrees;

    // Precompute the clusteringTypeMap outside of the loop to avoid redundant computations
    const std::unordered_map<std::string, int> clusteringTypeMap = {
        {"Complete", HCLUST_METHOD_COMPLETE},
        {"Average", HCLUST_METHOD_AVERAGE},
        {"Median", HCLUST_METHOD_MEDIAN},
        {"Single", HCLUST_METHOD_SINGLE} // Added "Single" to the map for consistency
    };
    std::string clusteringTypecurrentText = "Single";  // "Single", "Complete", "Average", "Median"
    int opt_method = clusteringTypeMap.at(clusteringTypecurrentText); 

    for (const auto& gene : returnGeneList) {
        QList<QStandardItem*> row;
        std::vector<float> numbers;
        std::map<QString, float> meanValuesForSpeciesMap;

        for (const auto& outerPair : map) {
            QString outerKey = outerPair.first;
            const std::map<QString, float>& innerMap = outerPair.second;
            auto it = innerMap.find(gene);
            float value = (it != innerMap.end()) ? it->second : 0.0f;
            numbers.push_back(value);
            meanValuesForSpeciesMap[outerKey] = value;
        }

        auto numOfLeaves = static_cast<int>(numOfSpecies);
        std::unique_ptr<double[]> distmat(condensedDistanceMatrix(numbers));

        std::unique_ptr<int[]> merge(new int[2 * (numOfLeaves - 1)]);
        std::unique_ptr<double[]> height(new double[numOfLeaves - 1]);
        hclust_fast(numOfLeaves, distmat.get(), opt_method, merge.get(), height.get());
        std::string newick = mergeToNewick(merge.get(), numOfLeaves);

        newickTrees.insert({ gene, {QString::fromStdString(newick), meanValuesForSpeciesMap} });

        row.push_back(new QStandardItem(gene)); //0
        row.push_back(new QStandardItem(""));  //1
        row.push_back(new QStandardItem(QString::number(-1.0)));  //2

        float meanV = calculateMean(numbers);
        // make meanV into 2 decimal places
        meanV = std::round(meanV * 100) / 100;
        row.push_back(new QStandardItem(QString::number(meanV)));//3

        auto it = geneCounter.find(gene);
        int count = (it != geneCounter.end()) ? it->second.size() : -1;
        row.push_back(new QStandardItem(QString::number(count)));//4

        QString speciesGeneAppearancesComb;
        if (it != geneCounter.end()) {
            for (const auto& str : it->second) {
                if (!speciesGeneAppearancesComb.isEmpty()) {
                    speciesGeneAppearancesComb += ";";
                }
                speciesGeneAppearancesComb += str;
            }
        }

        row.push_back(new QStandardItem(speciesGeneAppearancesComb));//5

        for (auto numb : numbers) {
            row.push_back(new QStandardItem(QString::number(numb)));
        }

        model->appendRow(row);
    }


    QString targetColor = "#fdb900";
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
                std::map <QString,float> speciesMeanMaps = pair.second.second;

                const char* string1 = targetNewick.c_str();
                const char* string2 = modifiedNewick.c_str();
                Tree t1;
                Tree t2;
                freopen("CON", "r", stdin);
                FILE* file1;
                FILE* file2;
                fopen_s(&file1, "file1.txt", "w");
                fputs(string1, file1);
                if (file1 != nullptr) {
                    fclose(file1);
                }

                fopen_s(&file2, "file2.txt", "w");
                fputs(string2, file2);
                if (file1 != nullptr) {
                    fclose(file2);
                }
                freopen("file1.txt", "r", stdin);
                t1.CreateTree();
                freopen("file2.txt", "r", stdin);
                t2.CreateTree();
                int sim = Calculate(&t1, &t2); 

                float similarity = 1.0 - static_cast<float>(sim) / static_cast<float>(numOfSpecies); 
                std::pair<QString, float>  temp;
                temp.first = createJsonTreeFromNewick(QString::fromStdString(modifiedNewick), leafnames, speciesMeanMaps);
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



    return QVariant::fromValue(model);

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


QString SettingsAction::createJsonTreeFromNewick(QString tree, std::vector<QString> leafnames, std::map <QString, float> speciesMeanValues)
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
                    meanValue = std::to_string(it->second);
                }
                else {
                    // Key not found, assign -1
                    meanValue = "-1";
                }

                jsonStream << "{\n\"color\": \"#000000\",\n\"hastrait\": true,\n\"iscollapsed\": false,\n\"branchLength\": 1.0,\n\"mean\": "<< meanValue <<", \n\"name\": \"" << species << "\"\n}";
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
    _startComputationTriggerAction.fromParentVariantMap(variantMap);
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
    _scatterplotEmbeddingColorOption.fromParentVariantMap(variantMap);
    _scatterplotEmbeddingPointsUMAPOption.fromParentVariantMap(variantMap);
    _selectedSpeciesVals.fromParentVariantMap(variantMap);
    _removeRowSelection.fromParentVariantMap(variantMap);
    _statusColorAction.fromParentVariantMap(variantMap);
    _typeofTopNGenes.fromParentVariantMap(variantMap);
    _usePreComputedTSNE.fromParentVariantMap(variantMap);
    _selectedCellClusterInfoBox.fromParentVariantMap(variantMap);
}

QVariantMap SettingsAction::toVariantMap() const
{
    QVariantMap variantMap = WidgetAction::toVariantMap();

    _geneNamesConnection.insertIntoVariantMap(variantMap);
    _startComputationTriggerAction.insertIntoVariantMap(variantMap);
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
    _scatterplotEmbeddingColorOption.insertIntoVariantMap(variantMap);
    _scatterplotEmbeddingPointsUMAPOption.insertIntoVariantMap(variantMap);
    _selectedSpeciesVals.insertIntoVariantMap(variantMap);
    _removeRowSelection.insertIntoVariantMap(variantMap);
    _statusColorAction.insertIntoVariantMap(variantMap);
    _typeofTopNGenes.insertIntoVariantMap(variantMap);
    _usePreComputedTSNE.insertIntoVariantMap(variantMap);
    _selectedCellClusterInfoBox.insertIntoVariantMap(variantMap);
    return variantMap;
}
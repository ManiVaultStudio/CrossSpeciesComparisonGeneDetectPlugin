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

    float sum = std::reduce(std::execution::par, v.begin(), v.end());
    float mean = sum / v.size();

    return mean;
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

    QStringList sortedList1 = list1;
    QStringList sortedList2 = list2;

    std::sort(sortedList1.begin(), sortedList1.end());
    std::sort(sortedList2.begin(), sortedList2.end());

    return sortedList1 == sortedList2;
}


int findIndex(const std::vector<std::seed_seq::result_type >& vec, int value) {
    auto it = std::find(vec.begin(), vec.end(), value);

    // If element was found
    if (it != vec.end()) {
        // Calculating the index
        int index = std::distance(vec.begin(), it);
        return index;
    }
    else {
        // Element not found
        return -1;
    }
}

SettingsAction::SettingsAction(CrossSpeciesComparisonGeneDetectPlugin& CrossSpeciesComparisonGeneDetectPlugin) :
    WidgetAction(&CrossSpeciesComparisonGeneDetectPlugin, "CrossSpeciesComparisonGeneDetectPlugin Settings"),
    _crossSpeciesComparisonGeneDetectPlugin(CrossSpeciesComparisonGeneDetectPlugin),
    _tableModel(this, "Table Model"),
    _selectedGene(this, "Selected Gene"),
    _filteringTreeDataset(this, "Filtering Tree Dataset"),
    _selectedRowIndex(this, "Selected Row Index"),
    _optionSelectionAction(*this),
    _startComputationTriggerAction(this, "Compute table"),
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
    _topNGenesFilter(this, "Top N Genes Filter", 10),
    _geneNamesConnection(this, "Gene Names Connection"),
    _createRowMultiSelectTree(this, "Create Row MultiSelect Tree"),
    _performGeneTableTsneAction(this, "Perform Gene Table TSNE"),
    _tsnePerplexity(this, "TSNE Perplexity"),
    _hiddenShowncolumns(this, "Hidden Shown Columns"),
    _scatterplotColorOption(this, "Scatterplot Color Option"),
    _selectedSpeciesVals(this, "Selected Species Vals")
{
    setSerializationName("CSCGDV:CrossSpeciesComparison Gene Detect Plugin Settings");
    _tableModel.setSerializationName("CSCGDV:Table Model");
    _selectedGene.setSerializationName("CSCGDV:Selected Gene");
    _mainPointsDataset  .setSerializationName("CSCGDV:Main Points Dataset");
    _embeddingDataset.setSerializationName("CSCGDV:Embedding Dataset");
    _speciesNamesDataset.setSerializationName("CSCGDV:Species Names Dataset");
    _clusterNamesDataset.setSerializationName("CSCGDV:Cluster Names Dataset");
    _filteredGeneNamesVariant.setSerializationName("CSCGDV:Filtered Gene Names");
    _topNGenesFilter.setSerializationName("CSCGDV:Top N Genes Filter");
    _filteringTreeDataset.setSerializationName("CSCGDV:Filtering Tree Dataset");
    _referenceTreeDataset.setSerializationName("CSCGDV:Reference Tree Dataset");
    _selectedRowIndex.setSerializationName("CSCGDV:Selected Row Index");
    _geneNamesConnection.setSerializationName("CSCGDV:Gene Names Connection");
    _selectedSpeciesVals.setSerializationName("CSCGDV:Selected Species Vals");
    _selectedGene.setDisabled(true);
    _selectedGene.setString("");
    _startComputationTriggerAction.setSerializationName("CSCGDV:Start Computation");
    _createRowMultiSelectTree.setSerializationName("CSCGDV:Create Row MultiSelect Tree");
    _performGeneTableTsneAction.setSerializationName("CSCGDV:Perform Gene Table TSNE");
    _tsnePerplexity.setSerializationName("CSCGDV:TSNE Perplexity");
    _tsnePerplexity.setMinimum(1);
    _tsnePerplexity.setMaximum(50);
    _tsnePerplexity.setValue(30);
    _hiddenShowncolumns.setSerializationName("CSCGDV:Hidden Shown Columns");
    _scatterplotColorOption.setSerializationName("CSCGDV:Scatterplot Color Option");
    _performGeneTableTsneAction.setChecked(false);
    _createRowMultiSelectTree.setDisabled(true);
    _selectedRowIndex.setDisabled(true);
    _selectedRowIndex.setString("");
    _scatterplotColorOption.initialize({"Species","Cluster","Expression"}, "Species");
    _filteringTreeDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
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

            _selectedSpeciesVals.setString("");
            auto pointsDataset = _mainPointsDataset.getCurrentDataset();
            auto embeddingDataset= _embeddingDataset.getCurrentDataset();
            auto speciesDataset = _speciesNamesDataset.getCurrentDataset();
            auto clusterDataset = _clusterNamesDataset.getCurrentDataset();
            auto referenceTreeDataset = _referenceTreeDataset.getCurrentDataset();
            auto filteringTreeDataset = _filteringTreeDataset.getCurrentDataset();

            if (_selectedPointsTSNEDataset.isValid())
            {
                _selectedPointsTSNEDataset->setSelectionIndices({});
            }

            bool isValid = false;
            QString datasetId = "";
            _geneNamesConnection.setString("");

            _clusterNameToGeneNameToExpressionValue.clear();
            //std::vector<QString> leafnames;
            if (pointsDataset.isValid() && speciesDataset.isValid() && clusterDataset.isValid() && referenceTreeDataset.isValid() && filteringTreeDataset.isValid() && embeddingDataset.isValid())
            {
                datasetId = referenceTreeDataset->getId();
                isValid = speciesDataset->getParent() == pointsDataset && clusterDataset->getParent()==pointsDataset && embeddingDataset->getParent() == pointsDataset;

                auto allSelectedIndices = pointsDataset->getSelectionIndices();
                auto rawEmbeddingDataset = mv::data().getDataset<Points>(embeddingDataset->getId());
                auto rawPointdata = mv::data().getDataset<Points>(pointsDataset->getId());
                auto allgeneList = rawPointdata->getDimensionNames();
                auto embeddingDimList = rawEmbeddingDataset->getDimensionNames();
                std::vector<int> embeddingGeneIndices;
                std::vector<int> allGeneIndices;
                for (int i = 0; i < allgeneList.size(); i++)
                {
                    allGeneIndices.push_back(i);
                }

                for (int i = 0; i < embeddingDimList.size(); i++)
                {
                    embeddingGeneIndices.push_back(i);
                }
                if (allSelectedIndices.size() > 0 && embeddingGeneIndices.size()>0)
                {


                    if (isValid)
                    {
                        auto speciesData = mv::data().getDataset<Clusters>(speciesDataset->getId());
                        auto clustersData = mv::data().getDataset<Clusters>(clusterDataset->getId());
                        auto clustersAll = clustersData->getClusters();
                        auto speciesAll = speciesData->getClusters();

                        std::map<QString, std::pair<QColor, std::vector<int>>> selctedClustersMap;
                        std::map<QString, std::pair<QColor, std::vector<int>>> selectedSpeciesMap;

                        if (!speciesAll.empty() && !clustersAll.empty())
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


                                
                                std::vector<float> resultContainerForSelectedPoints(allSelectedIndices.size()* allGeneIndices.size());
                                rawPointdata->populateDataForDimensions(resultContainerForSelectedPoints, allGeneIndices, allSelectedIndices);

                                QString datasetIdEmb = _selectedPointsDataset->getId();
                                int sizeofdatasetEmb = allSelectedIndices.size();
                                int dimofDatasetEmb = allGeneIndices.size();
                                populatePointData(datasetIdEmb, resultContainerForSelectedPoints, sizeofdatasetEmb, dimofDatasetEmb, allgeneList);




                                std::vector<float> resultContainerForSelectedEmbeddingPoints(allSelectedIndices.size()* embeddingGeneIndices.size());
                                rawEmbeddingDataset->populateDataForDimensions(resultContainerForSelectedEmbeddingPoints, embeddingGeneIndices, allSelectedIndices);

                                QString datasetId = _selectedPointsEmbeddingDataset->getId();
                                int sizeofdataset = allSelectedIndices.size();
                                int dimofDataset = embeddingGeneIndices.size();
                                populatePointData(datasetId, resultContainerForSelectedEmbeddingPoints, sizeofdataset, dimofDataset, embeddingDimList);
                                
                                
                                std::vector<float> resultContainerColorPoints(allSelectedIndices.size() * 1);
                                std::fill(resultContainerColorPoints.begin(), resultContainerColorPoints.end(), -1.0);

                                QString datasetIdExp = _tsneDatasetExpressionColors->getId();
                                int sizeofdatasetExp = allSelectedIndices.size();
                                int dimofDatasetExp = 1;
                                std::vector<QString> dimensionNamesExp = { "Expression" };

                                populatePointData(datasetIdExp, resultContainerColorPoints, sizeofdatasetExp, dimofDatasetExp, dimensionNamesExp);
                                
                                //_selectedPointsDataset->setData(resultContainerForSelectedPoints.data(), allSelectedIndices.size(), allGeneIndices.size());
                                //_selectedPointsDataset->setDimensionNames(allgeneList);
                                //events().notifyDatasetDataChanged(_selectedPointsDataset);

                                auto analysisPlugin = mv::plugins().requestPlugin<AnalysisPlugin>("tSNE Analysis", { _selectedPointsEmbeddingDataset });
                                if (!analysisPlugin) {
                                    qDebug() << "Could not find create TSNE Analysis";
                                    return;
                                }

                                int perplexity = std::min(static_cast<int>(allSelectedIndices.size()), _tsnePerplexity.getValue());
                                if (perplexity < 2)
                                {
                                    return;
                                }
                                if (perplexity != _tsnePerplexity.getValue())
                                {
                                    _tsnePerplexity.setValue(perplexity);
                                }

                                if (_selectedPointsTSNEDataset.isValid())
                                {
                                    auto runningAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDataset->findChildByPath("TSNE/TsneComputationAction/Running"));

                                    if (runningAction)
                                    {
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
                                _selectedPointsTSNEDataset = analysisPlugin->getOutputDataset();
                                if (_selectedPointsTSNEDataset.isValid())
                                {


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
                                                        


                                                        auto selectedColorType = _scatterplotColorOption.getCurrentText();
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

                                    auto startAction = dynamic_cast<TriggerAction*>(_selectedPointsTSNEDataset->findChildByPath("TSNE/TsneComputationAction/Start"));
                                    if (startAction) {

                                        startAction->trigger();

                                        analysisPlugin->getOutputDataset()->setSelectionIndices({});
                                    }

                                }


                            }
                            
                            for (auto& clusters : clustersAll)
                            {
                                auto clusterIndices = clusters.getIndices();
                                auto clusterName = clusters.getName();
                                auto clusterColor = clusters.getColor();
                                std::vector<int> filteredIndices; 
                                for (int i = 0; i < clusterIndices.size(); i++)
                                {

                                    int indexVal = findIndex(allSelectedIndices, clusterIndices[i]);
                                        if (indexVal != -1)
                                        {
                                            filteredIndices.push_back(indexVal);
                                        }

                                }
                                selctedClustersMap[clusterName] = { clusterColor, filteredIndices };
                            }
                            
                            for (auto& species : speciesAll)
                            {
                                auto speciesIndices = species.getIndices();
                                auto speciesName = species.getName();
                                auto speciesColor = species.getColor();

                                std::vector<int> filteredIndices;
                                for (int i = 0; i < speciesIndices.size(); i++)
                                {

                                    int indexVal = findIndex(allSelectedIndices, speciesIndices[i]);
                                    if (indexVal != -1)
                                    {
                                        filteredIndices.push_back(indexVal);
                                    }

                                }
                                selectedSpeciesMap[speciesName] = { speciesColor, filteredIndices };



                                //indices overlap between  speciesIndices and allSelectedIndices
                                std::vector<int> commonSelectedIndices;

                                std::sort(allSelectedIndices.begin(), allSelectedIndices.end());
                                std::sort(speciesIndices.begin(), speciesIndices.end());
                                std::set_intersection(allSelectedIndices.begin(), allSelectedIndices.end(), speciesIndices.begin(), speciesIndices.end(), std::back_inserter(commonSelectedIndices));

                                


                                for (int i = 0; i < allgeneList.size(); i++)
                                {
                                    auto geneName = allgeneList[i];
                                    auto geneIndex = { i };
                                    float meanValue = 0.0;
                                    if (commonSelectedIndices.size() > 0) {

                                        std::vector<float> resultContainerShort(commonSelectedIndices.size());
                                        rawPointdata->populateDataForDimensions(resultContainerShort, geneIndex, commonSelectedIndices);
                                        float shortMean = calculateMean(resultContainerShort);
                                        float fullMean = 0.0;
                                        if (_clusterGeneMeanExpressionMap[speciesName].find(geneName) == _clusterGeneMeanExpressionMap[speciesName].end())
                                        {
                                            std::vector<float> resultContainerFull(speciesIndices.size());
                                            rawPointdata->populateDataForDimensions(resultContainerFull, geneIndex, speciesIndices);
                                            fullMean = calculateMean(resultContainerFull);
                                            _clusterGeneMeanExpressionMap[speciesName][geneName] = fullMean;
                                        }
                                        else
                                        {
                                            fullMean = _clusterGeneMeanExpressionMap[speciesName][geneName];
                                        }
                                    
                                    if (fullMean != 0.0)
                                    {
                                        meanValue = shortMean/ fullMean;
                                    }
                                    //else
                                    //{
                                    //    meanValue = 0.0;
                                    //}
                                }

                                    _clusterNameToGeneNameToExpressionValue[speciesName][geneName]=meanValue;
                                }




                            }

                           

                            auto clusterColorDatasetId= _tsneDatasetClusterColors->getId();
                            auto speciesColorDatasetId = _tsneDatasetSpeciesColors->getId();
                            
                            populateClusterData(speciesColorDatasetId, selectedSpeciesMap);
                            populateClusterData(clusterColorDatasetId, selctedClustersMap);

                            QVariant geneListTable = findTopNGenesPerCluster(_clusterNameToGeneNameToExpressionValue, _topNGenesFilter.getValue(), datasetId, 1.0);

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
                    }

                    else
                    {
                        qDebug() << "Invalid child datasets";
                    }


                }
                else
                {
                    qDebug() << "No points selected or no dimensions present";
                }
            }

            else
            {
                qDebug() << "Invalid datasets";
            }


        };
       
    connect(&_startComputationTriggerAction, &TriggerAction::triggered, this, updateGeneFilteringTrigger);
    const auto updateCreateRowMultiSelectTreeTrigger = [this]() -> void {
        
        if (_filteringTreeDataset.getCurrentDataset().isValid())
        {
            auto treeDataset = mv::data().getDataset<CrossSpeciesComparisonTree>(_filteringTreeDataset.getCurrentDataset().getDatasetId());

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
        else
        {
            qDebug() << "Tree dataset is not valid";
        }

        };

    connect(&_createRowMultiSelectTree, &TriggerAction::triggered, this, updateCreateRowMultiSelectTreeTrigger);

    const auto updateMainPointsDataset = [this]() -> void {

        _clusterGeneMeanExpressionMap.clear();
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

    const auto updateScatterplotColor = [this]() -> void {
        auto selectedColorType= _scatterplotColorOption.getCurrentText();
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
                                



                                auto selectedColorType = _scatterplotColorOption.getCurrentText();
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
    connect(&_scatterplotColorOption, &OptionAction::currentIndexChanged, this, updateScatterplotColor);
    
    
    const auto updateEmbeddingDataset = [this]() -> void {


        };
    connect(&_embeddingDataset, &DatasetPickerAction::currentIndexChanged, this, updateEmbeddingDataset);  
    
    

}
QVariant SettingsAction::findTopNGenesPerCluster(const std::map<QString, std::map<QString, float>>& map, int n, QString datasetId, float treeSimilarityScore) {

    if (map.empty() || n <= 0) {
        return QVariant();
    }

    QSet<QString> geneList;
    QStringList returnGeneList;
    std::map<QString, std::vector<QString>> geneAppearanceCounter;

    for (const auto& outerPair : map) {
        // Convert map to vector of pairs
        std::vector<std::pair<QString, float>> geneExpressionVec(outerPair.second.begin(), outerPair.second.end());

        // Sort the vector in descending order based on the expression value
        std::sort(geneExpressionVec.begin(), geneExpressionVec.end(), [](const auto& a, const auto& b) {
            return a.second > b.second;
            });

        // Add the top n genes to the geneList and initialize their count in geneAppearanceCounter
        for (int i = 0; i < std::min(n, static_cast<int>(geneExpressionVec.size())); ++i) {
            geneList.insert(geneExpressionVec[i].first);
        }
    }

    // Convert QSet<QString> geneList to QStringList returnGeneList
    returnGeneList = QStringList(geneList.begin(), geneList.end());

    // Increment count for each gene in geneList for the geneAppearanceCounter for top n genes
    for (const auto& outerPair : map) {
        // Convert map to vector of pairs
        std::vector<std::pair<QString, float>> geneExpressionVec(outerPair.second.begin(), outerPair.second.end());

        // Sort the vector in descending order based on the expression value
        std::sort(geneExpressionVec.begin(), geneExpressionVec.end(), [](const auto& a, const auto& b) {
            return a.second > b.second;
            });

        for (int i = 0; i < std::min(n, static_cast<int>(geneExpressionVec.size())); ++i) {
            // If geneExpressionVec[i].first is present in the key of geneAppearanceCounter, push the outerpair.first in the vector
            if (geneList.contains(geneExpressionVec[i].first)) {
                geneAppearanceCounter[geneExpressionVec[i].first].push_back(outerPair.first);
            }

        }
    }

    //qDebug() << "***********Insert location\n";
    //for (auto& pair : geneAppearanceCounter) {
    //    std::cout << "Gene: " << pair.first.toStdString() << ", Count: " << pair.second << std::endl;
    //}

    QVariant returnValue = createModelFromData(returnGeneList, map, datasetId, treeSimilarityScore, geneAppearanceCounter, n);

    return returnValue;
}

QVariant SettingsAction::createModelFromData(const QStringList& returnGeneList, const std::map<QString, std::map<QString, float>>& map, const QString& treeDatasetId, const float& treeSimilarityScore, const std::map<QString, std::vector<QString>>& geneCounter, const int& n) {

    if (returnGeneList.isEmpty() || map.empty()) {
        return QVariant();
    }

    QStandardItemModel* model = new QStandardItemModel();
    int numOfSpecies = map.size();
    QStringList initColumnNames = { "ID", "Newick tree","Tree Similarity with Reference Tree", "Total Mean","Top Gene Appearances/" + QString::number(numOfSpecies) + " Species", "Gene Apearance Species Names" };
    model->setHorizontalHeaderLabels(initColumnNames);
    
    std::map<QString, std::map<QString, float>>::const_iterator it = map.begin();
    for (int i = 0 + initColumnNames.size(); i < numOfSpecies + initColumnNames.size(); i++, it++) {
        QString headerTitle = it->first;
        //headerTitle.replace("_", " ");
        //headerTitle = QString("Mean ") + headerTitle;
        model->setHorizontalHeaderItem(i, new QStandardItem(headerTitle));
    }
    QStringList headers;
    for (int i = 0; i < model->columnCount(); ++i) {
        headers.push_back(model->horizontalHeaderItem(i)->text());
    }
    _hiddenShowncolumns.setOptions(headers);

    QStringList selectedHeaders = { headers[0], headers[2], headers[3], headers[4] };// , headers[5]};
    _hiddenShowncolumns.setSelectedOptions(selectedHeaders);

    std::map<QString, std::pair<QString, std::map<QString, float>>> newickTrees;
    //std::map<QString, std::map<QString, float>> meanExpressionValues;
    for (auto gene : returnGeneList)
    {
        QList<QStandardItem*> row;
        std::vector<float> numbers;
        std::map<QString, float> meanValuesForSpeciesMap;
        //[speciesName][geneName]=meanValue;
        for (const auto& outerPair : map) {
            QString outerKey = outerPair.first;
            const std::map<QString, float>& innerMap = outerPair.second;
            try {
                numbers.push_back(innerMap.at(gene));
                meanValuesForSpeciesMap[outerKey] = innerMap.at(gene);
            }
            catch (const std::out_of_range& e) {
                numbers.push_back(0);
                meanValuesForSpeciesMap[outerKey] = 0;
            }



        }

        const std::unordered_map<std::string, int> clusteringTypeMap = {
    {"Complete", HCLUST_METHOD_COMPLETE},
    {"Average", HCLUST_METHOD_AVERAGE},
    {"Median", HCLUST_METHOD_MEDIAN}
        };
        std::string clusteringTypecurrentText = "Single";  //"Single","Complete", "Average","Median"
        int opt_method = clusteringTypeMap.count(clusteringTypecurrentText) ? clusteringTypeMap.at(clusteringTypecurrentText) : HCLUST_METHOD_SINGLE;

        auto numOfLeaves = numOfSpecies;
        double* distmat = new double[(numOfLeaves * (numOfLeaves - 1)) / 2];


        distmat = condensedDistanceMatrix(numbers);

        int* merge = new int[2 * (numOfLeaves - 1)];
        double* height = new double[numOfLeaves - 1];
        hclust_fast(numOfLeaves, distmat, opt_method, merge, height);
        std::string newick = mergeToNewick(merge, numOfLeaves);
        int totalChars = newick.length();
        //std::cout << "\nOriginal Newick format: " << newick << std::endl;
        //add gene and newick to newickTrees
       // newickTrees.insert(std::make_pair(gene, QString::fromStdString(newick)));
        newickTrees.insert(std::make_pair(gene, std::make_pair(QString::fromStdString(newick), meanValuesForSpeciesMap)));

        delete[] distmat;
        delete[] merge;
        delete[] height;

        //Statistics stats = calculateStatistics(numbers);

        row.push_back(new QStandardItem(gene));
        //row.push_back(new QStandardItem(QString::number(stats.variance)));
        row.push_back(new QStandardItem(""));
        //row.push_back(new QStandardItem(QString::number(-1)));
        row.push_back(new QStandardItem()), row.back()->setData(-1, Qt::DisplayRole), row.back()->setData(-1, Qt::UserRole);
        float meanV= calculateMean(numbers);
        row.push_back(new QStandardItem()), row.back()->setData(meanV, Qt::DisplayRole), row.back()->setData(meanV, Qt::UserRole);
        QString key = gene;
        //qDebug() << "\n**Trying to find key:" << gene << "\n";
        int count = -1;
        auto it = geneCounter.find(key);
        if (it != geneCounter.end()) {
            //qDebug()<< "Species counter"<< key << "found.\n";
            //qDebug()<< "it->second"<< it->second << "found.\n";
            //qDebug()<< "(it->second).size()"<< (it->second).size() << "found.\n";
            count = (it->second).size();
            //row.push_back(new QStandardItem(QString::number(count)));
            row.push_back(new QStandardItem()), row.back()->setData(count, Qt::DisplayRole), row.back()->setData(count, Qt::UserRole);
        }
        else {
            qDebug() << "Key " << gene << "not found.\n";
            //row.push_back(new QStandardItem(QString::number(-1)));
            row.push_back(new QStandardItem()), row.back()->setData(count, Qt::DisplayRole), row.back()->setData(count, Qt::UserRole);
        }

        //row.push_back(new QStandardItem(QString::number(stats.stdDeviation)));
        //row.push_back(new QStandardItem(QString::number(stats.mean)));
        QString speciesGeneAppearancesComb;
        for (const auto& str : it->second) {
            if (!speciesGeneAppearancesComb.isEmpty()) {
                speciesGeneAppearancesComb += ";";
            }
            speciesGeneAppearancesComb += str;
        }

        //row.push_back(new QStandardItem(speciesGeneAppearancesComb));
        row.push_back(new QStandardItem()), row.back()->setData(speciesGeneAppearancesComb, Qt::DisplayRole), row.back()->setData(count, Qt::UserRole);
        for (auto numb : numbers)
        {
            //row.push_back(new QStandardItem(QString::number(numb)));
            row.push_back(new QStandardItem()), row.back()->setData(numb, Qt::DisplayRole), row.back()->setData(numb, Qt::UserRole);
        }

        // Create a new item for the vertical header
        //QStandardItem* verticalHeaderItem = new QStandardItem(QString::number(geneCounter[gene]));

        // Add the new item to the model's vertical header
       // model->setVerticalHeaderItem(model->rowCount(), verticalHeaderItem);

        // Add the row to the model
        model->appendRow(row);
    }

    //qDebug() << "***********Access location\n";
    //for (auto& pair : geneCounter) {
    //    std::cout << "Gene: " << pair.first.toStdString() << ", Count: " << pair.second << std::endl;
    //}

    //print newickTrees
    //for (auto& pair : newickTrees) {
    //    std::cout << "Gene: " << pair.first.toStdString() << ", Newick: " << pair.second.toStdString() << std::endl;
    //}



    //check which newick trees are exactly similar to "(((((((20,(((24,((25,9),(12,11))),(19,15)),(22,21))),((1,17),23)),((2,18),(8,6))),(14,10)),(3,16)),(7,5)),(13,4));" and add color "#00a2ed" to the genes
    // Define the target newick tree and color
    //std::string targetNewick = "((((((20,(((((18,24),(19,17)),16),(11,25)),(9,6))),(23,4)),((3,1),(15,(22,((13,2),7))))),(10,5)),(8,21)),(12,14));";
    QString targetColor = "#fdb900";
    std::string targetNewick = "";
    QStringList fullTreeNames;

    //iterate std::map<QString, std::map<QString, float>> map and fill keys in std::vector<QString> leafNames 
    std::vector<QString> leafnames;
    for (const auto& outerPair : map) {
        leafnames.push_back(outerPair.first);
    }



    std::map<QString, std::pair<QString, float>> treeSimilarities;
    if (treeDatasetId != "")
    {
        auto fullTreeData = mv::data().getDataset<CrossSpeciesComparisonTree>(treeDatasetId);
        if (fullTreeData.isValid())
        {
            fullTreeNames = fullTreeData->getTreeLeafNames();
            auto treeData = fullTreeData->getTreeData();

            QJsonDocument jsonDoc(treeData);
            auto temp = jsonDoc.toJson(QJsonDocument::Compact).toStdString();
            auto jsonTree = nlohmann::json::parse(temp);
            targetNewick = jsonToNewick(jsonTree, leafnames);
            targetNewick += ";";  // End of Newick string
        }

    }

    /*

    auto jsonTree = nlohmann::json::parse(jsonString);
    std::string targetNewick = jsonToNewick(jsonTree, leafnames);
    targetNewick += ";";  // End of Newick string

    */
    if (fullTreeNames.size() > 0 && leafnames.size() > 0 && targetNewick != "")
    {

        //convert  std::vector<QString> to QStringList leafnames
        QStringList copyleafNames;
        for (auto& leaf : leafnames) {
            copyleafNames.push_back(leaf);
        }



        if (areSameIgnoreOrder(fullTreeNames, copyleafNames)) {
            // Iterate over the newickTrees map
            for (auto& pair : newickTrees) {
                /*
                        qDebug() << "\n*****************";
                        qDebug() << "First tree: " << pair.second;
                        qDebug() << "Second tree: " << QString::fromStdString(targetNewick);
                        qDebug() << "*****************\n";
                        */
                        //add a ";" to the end of the string pair.second.toStdString()
                std::string modifiedNewick = pair.second.first.toStdString();
                std::map <QString,float> speciesMeanMaps = pair.second.second;

                const char* string1 = targetNewick.c_str();
                const char* string2 = modifiedNewick.c_str();

                //const char* string1 = "(((((((20,(((24,((25,9),(12,11))),(19,15)),(22,21))),(23,(1,17))),((2,18),(8,6))),(14,10)),(3,16)),(7,5)),(13,4));";
                //const char* string2 = "(((((((20,(((24,((25,9),(12,11))),(19,15)),(22,21))),((1,17),23)),((2,18),(8,6))),(14,10)),(16,3)),(7,5)),(13,4));";

            // Create two Tree objects
                Tree t1;
                Tree t2;

                // Change the standard input to read from the strings
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

                // Read tree structures from the strings
                freopen("file1.txt", "r", stdin);
                t1.CreateTree();
                freopen("file2.txt", "r", stdin);
                t2.CreateTree();

                // Calculate and print the similarity
                int sim = Calculate(&t1, &t2); // sim is the minimum number of leaves to remove to make the trees isomorphic

                /* qDebug() << "\n*****************\n"
                     << "First tree: " << string1
                     << "\nSecond tree: " << string2
                     << "\nSimvalue: " << sim
                     << "\n*****************\n"; */

                     //qDebug()<<"\n****Simvalue: "<<sim<<"****\n";

                     // If the current newick tree is the same as the target

                float similarity = 1.0 - static_cast<float>(sim) / static_cast<float>(numOfSpecies); //the similarity between two Newick trees,


                //insert pair.first modifiedNewick similarity to treeSimilarities
                std::pair<QString, float>  temp;
                temp.first = createJsonTreeFromNewick(QString::fromStdString(modifiedNewick), leafnames, speciesMeanMaps);
                temp.second = similarity;
                treeSimilarities.insert(std::make_pair(pair.first, temp));
                /*
                 1.	Calculate(&t1, &t2) is a function that takes two trees in Newick format and returns the minimum number of leaves that need to be removed to make them isomorphic.
        2.	sim is the result of this calculation.
        3.	numOfSpeciesLeaves is presumably the total number of leaves in the tree (or in both trees if they have the same number of leaves).
        4.	static_cast<float>(sim) / static_cast<float>(numOfSpeciesLeaves) calculates the proportion of leaves that need to be removed to make the trees isomorphic.
        5.	1.0 - static_cast<float>(sim) / static_cast<float>(numOfSpeciesLeaves) then subtracts this proportion from 1 to give the proportion of leaves that do not need to be removed, which can be interpreted as a measure of similarity between the trees.
                 if no leaves need to be removed (i.e., the trees are already isomorphic), sim will be 0, and the similarity will be 1.0. If all leaves need to be removed, sim will be equal to numOfSpeciesLeaves, and the similarity will be 0.
                */
                /*
                int x= (1-treeSimilarityScore)*numOfSpecies;

                if (sim <= x) {
                    // Find the corresponding gene in the model
                    QList<QStandardItem*> items = model->findItems(pair.first);
                    for (auto& item : items) {
                        // Get the row of the item
                        int row = item->row();
                        // Set the background color of each item in the row
                        for (int i = 0; i < model->columnCount(); i++) {
                            QStandardItem* itemInRow = model->item(row, i);
                            if (itemInRow) {
                                itemInRow->setBackground(QBrush(QColor(targetColor)));
                            }
                        }
                    }
                }

                */
            }

        }
    }

    //check which newick trees are the same and group them together
    /*
    std::map<QString, std::vector<QString>> clusteringMap;

    for (auto it = newickTrees.begin(); it != newickTrees.end(); ++it) {
        QString currentNewick = it->second;
        QString currentGene = it->first;
        if (clusteringMap.empty()) {
            clusteringMap[currentNewick] = { currentGene };
        }
        else {
            bool found = false;
            for (auto& cluster : clusteringMap) {
                QString clusterNewick = cluster.first;
                if (currentNewick == clusterNewick) {
                    cluster.second.push_back(currentGene);
                    found = true;
                    break;
                }
            }
            if (!found) {
                clusteringMap[currentNewick] = { currentGene };
            }
        }
    }
    */
    //QStringList colorCodes = {"#8dd3c7" ,"#ffffb3"}

    //colors
    // left right color codes https://encycolorpedia.com/00a2ed #ff5d12 and #00a2ed  or #0fb1e0 and #f04e1f
    //testing
    //clusteringMap[clusteringMap.begin()->first] = { "VSTM5", "REC8", "TPGS2" }; 
    //clusteringMap[(++clusteringMap.begin())->first] = { "FGD6", "ANGEL1","FANCC" };


    //print clusteringMap    '
    /*for (auto& cluster : clusteringMap) {
        QString newick = cluster.first;
        std::vector<QString> genes = cluster.second;

        if (genes.size() > 1) {
            std::cout << "Newick: " << newick.toStdString() << std::endl;
            std::cout << "Genes: ";
            for (auto& gene : genes) {
                std::cout << gene.toStdString() << ", ";
            }
            std::cout << std::endl;
        }

    }
    */
    /*
    QStringList colorCodes = { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928" };
    int colorIndex = 0;
    for (auto& cluster : clusteringMap) {
        QString newick = cluster.first;
        std::vector<QString> genes = cluster.second;

        if (genes.size() > 1) {
            for (auto& gene : genes) {
                for (int i = 0; i < model->rowCount(); i++) {
                    if (model->item(i, 0)->text() == gene) {
                        for (int j = 0; j < model->columnCount(); j++) {
                            model->item(i, j)->setBackground(QBrush(QColor(colorCodes[colorIndex])));
                        }
                    }
                }
            }
            colorIndex = (colorIndex + 1) % colorCodes.size();
        }
    }

    */

    //based on first column string value from  model, update the 4th column vaLUE   from treeSimilarities
    for (int i = 0; i < model->rowCount(); i++) {
        QString gene = model->item(i, 0)->text();
        //qDebug() <<"Gene: " << gene;
        //qDebug() << "Tree Similarity: " << treeSimilarities[gene];
        auto it = treeSimilarities.find(gene);
        auto map = it->second;
        auto similarity = map.second;
        auto newick = map.first;

        if (it != treeSimilarities.end()) {

            model->item(i, 1)->setText(newick);
            //model->item(i, 2)->setText(QString::number(similarity));

            model->item(i, 2)->setData(similarity, Qt::DisplayRole);
            model->item(i, 2)->setData(similarity, Qt::UserRole);


        }
    }



    return QVariant::fromValue(model);

}
void SettingsAction::populatePointData(QString& datasetId, std::vector<float>& pointVector, int& numPoints, int& numDimensions, std::vector<QString>& dimensionNames)
{
    auto pointDataset = mv::data().getDataset<Points>(datasetId);
    if (pointDataset.isValid())
    {
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

double* SettingsAction::condensedDistanceMatrix(std::vector<float>& items) {
    int n = items.size();
    double* distmat = new double[(n * (n - 1)) / 2];
    int k = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            distmat[k] = std::abs(items[i] - items[j]);
            ++k;
        }
    }

    /*std::cout << "Distance matrix: " << std::endl;
    int index = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            std::cout << "Distance " << i << " value: " << items[i] << " and " << j << " value: " << items[j] << ": " << distmat[index++] << std::endl;
        }
    }*/
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
    _filteringTreeDataset.fromParentVariantMap(variantMap);
    _referenceTreeDataset.fromParentVariantMap(variantMap);
    _selectedRowIndex.fromParentVariantMap(variantMap);
    _performGeneTableTsneAction.fromParentVariantMap(variantMap);
    _tsnePerplexity.fromParentVariantMap(variantMap);
    _hiddenShowncolumns.fromParentVariantMap(variantMap);
    _scatterplotColorOption.fromParentVariantMap(variantMap);
    _selectedSpeciesVals.fromParentVariantMap(variantMap);
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
    _filteringTreeDataset.insertIntoVariantMap(variantMap);
    _referenceTreeDataset.insertIntoVariantMap(variantMap);
    _selectedRowIndex.insertIntoVariantMap(variantMap);
    _performGeneTableTsneAction.insertIntoVariantMap(variantMap);
    _tsnePerplexity.insertIntoVariantMap(variantMap);
    _hiddenShowncolumns.insertIntoVariantMap(variantMap);
    _scatterplotColorOption.insertIntoVariantMap(variantMap);
    _selectedSpeciesVals.insertIntoVariantMap(variantMap);
    return variantMap;
}
#include "CrossSpeciesComparisonGeneDetectPlugin.h"

#include <event/Event.h>
#include <CrossSpeciesComparisonTreeData.h>
#include <DatasetsMimeData.h>
#include <QHeaderView> 
#include <QDebug>
#include <QMimeData>
#include <QShortcut>
#include <QSplitter>
#include <QRandomGenerator>
#include <QColor>
#include <QJsonArray> 
#include <unordered_set>
#include <cmath>
#include <algorithm>
//#include<QTooltip>
#include <QRegularExpression> 
#ifdef _WIN32
#include <execution>
#endif
Q_PLUGIN_METADATA(IID "studio.manivault.CrossSpeciesComparisonGeneDetectPlugin")

using namespace mv;


void applyLogTransformation(std::vector<float>& values) {
    

#ifdef _WIN32
    std::transform(std::execution::par, values.begin(), values.end(), values.begin(),
        [](float value) { return std::log(value + 1); });
#else
    std::transform(values.begin(), values.end(), values.begin(),
        [](float value) { return std::log(value + 1); });
#endif

}

void clearTableSelection(QTableView* tableView) {
    if (tableView && tableView->selectionModel()) {
        // Clear the current selection
        tableView->clearSelection();

        // Temporarily disable the selection mode to remove highlight
        QAbstractItemView::SelectionMode oldMode = tableView->selectionMode();
        tableView->setSelectionMode(QAbstractItemView::NoSelection);

        // Clear the current index
        tableView->selectionModel()->setCurrentIndex(QModelIndex(), QItemSelectionModel::NoUpdate);

        // Restore the original selection mode
        tableView->setSelectionMode(oldMode);

        // Update the view to ensure changes are reflected
        tableView->update();
    }
    else {
        qDebug() << "TableView or its selection model is null";
    }
}
void updateRowVisibility(const QSet<QString>& uniqueReturnGeneList, QTableView* geneTableView, QSortFilterProxyModel* proxyModel) {

    for (int i = 0; i < proxyModel->rowCount(); i++) {
        QString geneName = proxyModel->index(i, 0).data().toString();
        int geneAppearanceSpeciesNamesCount = proxyModel->index(i, 1).data().toInt();

        if (!uniqueReturnGeneList.contains(geneName) || geneAppearanceSpeciesNamesCount < 1) {
            geneTableView->hideRow(i);
        }
        else {
            geneTableView->showRow(i);
        }
    }
}

void makeAllRowsVisible(QTableView* geneTableView, QSortFilterProxyModel* proxyModel) {
    int rowCount = proxyModel->rowCount();
    for (int row = 0; row < rowCount; ++row) {
        geneTableView->showRow(row);
    }
}


std::map<QString, SpeciesDetailsStats> convertToStatisticsMap(const QString& formattedStatistics) {
    std::map<QString, SpeciesDetailsStats> statisticsMap;

    QStringList speciesStatsList = formattedStatistics.split(";", Qt::SkipEmptyParts); // Qt 5.14 and later

    // qDebug() << speciesStatsList;
    QRegularExpression regex(R"(Species: (.*), Rank: (\d+), AbundanceTop: (\d+),  AbundanceMiddle: (\d+),  CountAbundanceNumerator: (\d+), MeanSelected: ([\d.]+), CountSelected: (\d+), MeanNotSelected: ([\d.]+), CountNotSelected: (\d+))");

    for (const QString& speciesStats : speciesStatsList) {
        QRegularExpressionMatch match = regex.match(speciesStats.trimmed());
        if (match.hasMatch()) {
            QString species = match.captured(1);
            SpeciesDetailsStats stats = {
                match.captured(2).toInt(),  // Rank
                match.captured(3).toInt(),  // AbundanceTop
                match.captured(4).toInt(),  // AbundanceMiddle
                match.captured(5).toInt(),  // CountAbundanceNumerator
                match.captured(6).toFloat(),  // MeanSelected
                match.captured(7).toInt(),  // CountSelected
                match.captured(8).toFloat(),  // MeanNotSelected
                match.captured(9).toInt()  // CountNotSelected
            };
            statisticsMap[species] = stats;
        }
    }

    return statisticsMap;
}



CrossSpeciesComparisonGeneDetectPlugin::CrossSpeciesComparisonGeneDetectPlugin(const PluginFactory* factory) :
    ViewPlugin(factory),
    _settingsAction(*this)
{

}

void CrossSpeciesComparisonGeneDetectPlugin::init()
{

    const auto updateSelectedRowIndex = [this]() -> void
        {

            if (_settingsAction.getFilteringEditTreeDatasetAction().getCurrentDataset().isValid())
            {
                auto treeDataset = mv::data().getDataset<CrossSpeciesComparisonTree>(_settingsAction.getFilteringEditTreeDatasetAction().getCurrentDataset().getDatasetId());
              
                QStringList selectedRowsStrList = _settingsAction.getSelectedRowIndexAction().getString().split(",");
                QList<int> selectedRows;
                for (const QString& str : selectedRowsStrList) {
                    selectedRows << str.toInt();
                }

                if (selectedRows.size()==1)
                {
                    int selectedRow = selectedRows[0];
                    if (treeDataset.isValid() && _settingsAction.getGeneTableView() && selectedRow >= 0)
                    {
                        QString treeData = _settingsAction.getGeneTableView()->model()->index(selectedRow, 2).data().toString();
                        //qDebug()<< "Tree data: " << treeData;
                        if (!treeData.isEmpty())
                        {

                            QJsonObject valueStringReference = QJsonDocument::fromJson(treeData.toUtf8()).object();
                            if (!valueStringReference.isEmpty())
                            {
                                treeDataset->setTreeData(valueStringReference);
                                events().notifyDatasetDataChanged(treeDataset);
                                //QString firstColumnValue = _settingsAction.getGeneTableView()->model()->index(selectedRow, 0).data().toString();
                               // _settingsAction.getGeneNamesConnection().setString(firstColumnValue);
                            }
                        }
                    }
                }
                if (selectedRows.size() > 1)
                {
                    _settingsAction.getCreateRowMultiSelectTree().setEnabled(true);
                }
                else
                {
                    _settingsAction.getCreateRowMultiSelectTree().setDisabled(true);
                }
                QStringList firstColumnValues;
                for (int row : selectedRows) {
                    firstColumnValues << _settingsAction.getGeneTableView()->model()->index(row, 0).data().toString();
                }
                QString firstColumnValue = firstColumnValues.join("*%$@*@$%*");
                //_settingsAction.getGeneNamesConnection().setString(firstColumnValue);


            }

        };

    connect(&_settingsAction.getSelectedRowIndexAction(), &StringAction::stringChanged, this, updateSelectedRowIndex);

    const auto updateSelectedGene = [this]() -> void
        {


        };

    connect(&_settingsAction.getSelectedGeneAction(), &StringAction::stringChanged, this, updateSelectedGene);

    const auto removeRowSelectionTable = [this]() -> void
        {
            auto statusString = _settingsAction.getStatusColorAction().getString();
            if (_settingsAction.getGeneTableView() && _settingsAction.getGeneTableView()->selectionModel()) {
                // Clear the current index if there's no selection
                _settingsAction.getGeneTableView()->clearSelection();

                // Temporarily disable the selection mode to remove highlight
                QAbstractItemView::SelectionMode oldMode = _settingsAction.getGeneTableView()->selectionMode();
                _settingsAction.getGeneTableView()->setSelectionMode(QAbstractItemView::NoSelection);

                // Clear the current index
                _settingsAction.getGeneTableView()->selectionModel()->setCurrentIndex(QModelIndex(), QItemSelectionModel::NoUpdate);

                // Restore the original selection mode
                _settingsAction.getGeneTableView()->setSelectionMode(oldMode);
                // Update the view to ensure changes are reflected
                _settingsAction.getGeneTableView()->update();
                _settingsAction.getSelctedSpeciesVals().setString("");


                if (_settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset().isValid() && _settingsAction.getBottomClusterNamesDataset().getCurrentDataset().isValid())
                {

                    auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
                    mv::gui::DatasetPickerAction* colorDatasetPickerAction;
                    mv::gui::DatasetPickerAction* pointDatasetPickerAction;
                    mv::gui::ViewPluginSamplerAction* samplerActionAction;


                    if (scatterplotViewFactory) {
                        for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                            if (plugin->getGuiName() == "Scatterplot Embedding View") {
                                pointDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Position"));
                                if (pointDatasetPickerAction) {
                                    pointDatasetPickerAction->setCurrentText("");

                                    pointDatasetPickerAction->setCurrentDataset(_settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset());

                                    colorDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Color"));
                                    if (colorDatasetPickerAction)
                                    {
                                        colorDatasetPickerAction->setCurrentText("");
                                        colorDatasetPickerAction->setCurrentDataset(_settingsAction.getBottomClusterNamesDataset().getCurrentDataset());

                                    }
                                    samplerActionAction = plugin->findChildByPath<mv::gui::ViewPluginSamplerAction>("Sampler");

                                    if (samplerActionAction)
                                    {
                                        samplerActionAction->setTooltipGeneratorFunction([this](const ViewPluginSamplerAction::SampleContext& toolTipContext) -> QString {
                                            QString clusterDatasetId = _settingsAction.getSpeciesNamesDataset().getCurrentDataset().getDatasetId();
                                            return _settingsAction.generateTooltip(toolTipContext, clusterDatasetId,true, "GlobalPointIndices");
                                            });
                                    }

                                }
                            }
                        }
                    }
                    _settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset()->setSelectionIndices(_settingsAction.getSelectedIndicesFromStorage());
                    mv::events().notifyDatasetDataSelectionChanged(_settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset());
                }
            }
            else {
                qDebug() << "TableView or its selection model is null";
            }



            if ( _settingsAction.getSelectedPointsTSNEDatasetForGeneTable().isValid())
            {

                {
                    _settingsAction.getSelectedPointsTSNEDatasetForGeneTable()->setSelectionIndices({});
                    mv::events().notifyDatasetDataSelectionChanged(_settingsAction.getSelectedPointsTSNEDatasetForGeneTable());

                }
            }




            _settingsAction.getRemoveRowSelection().setDisabled(true);
            //_settingsAction.enableDisableButtonsAutomatically();
            _settingsAction.getStatusColorAction().setString(statusString);
            selectedCellStatisticsStatusBarRemove();
            selectedCellCountStatusBarAdd();
            _settingsAction.getSpeciesExplorerInMap().setSelectedOptions({});
            _settingsAction.getClearRightClickedCluster().trigger();
        };

    connect(&_settingsAction.getRemoveRowSelection(), &TriggerAction::triggered, this, removeRowSelectionTable);


    const auto updateSpeciesExplorerInMap = [this]() -> void
        {
          
            geneExplorer();
            _settingsAction.enableDisableButtonsAutomatically();
            clearTableSelection(_settingsAction.getSelectionDetailsTable());
        };

    connect(&_settingsAction.getSpeciesExplorerInMapTrigger(), &TriggerAction::triggered, this, updateSpeciesExplorerInMap);


    const auto updateRevertRowSelectionChangesToInitial = [this]() -> void
        {
            // Check if _settingsAction is valid and initialized
            if (&_settingsAction) {
                QString selectedSpeciesVals = _settingsAction.getSelctedSpeciesVals().getString();

                // Check if the string is not empty
                if (!selectedSpeciesVals.isEmpty()) {
                    QStringList autoSpecies = selectedSpeciesVals.split(" @%$,$%@ ");

                    // Further check if autoSpecies is not just a list with an empty string
                    if (!(autoSpecies.size() == 1 && autoSpecies.first().isEmpty())) {
                        if (_settingsAction.getSpeciesExplorerInMap().getNumberOfOptions() > 0)
                        {
                            _settingsAction.getSpeciesExplorerInMap().setSelectedOptions(autoSpecies);
                            //_settingsAction.getRevertRowSelectionChangesToInitial().setDisabled(true);
                            geneExplorer();
                        }
                        
                    }

                }

            }
            //_settingsAction.getRevertRowSelectionChangesToInitial().setDisabled(true);
            clearTableSelection(_settingsAction.getSelectionDetailsTable());
            _settingsAction.getClearRightClickedCluster().trigger();
        };

    connect(&_settingsAction.getRevertRowSelectionChangesToInitial(), &TriggerAction::triggered, this, updateRevertRowSelectionChangesToInitial);


    const auto updateSelctedSpeciesVals = [this]() -> void
        {
            _settingsAction.enableDisableButtonsAutomatically();

        };

    connect(&_settingsAction.getSelctedSpeciesVals(), &StringAction::stringChanged, this, updateSelctedSpeciesVals);


    const auto updateTableModel = [this]() -> void
        {
            //_settingsAction.setErrorOutFlag(false);
            modifyListData();
            //if (_settingsAction.getErroredOutFlag())
            //{
                //_settingsAction.getStatusColorAction().setString("E");
                
            //}
            //else
            //{
                _settingsAction.getRemoveRowSelection().trigger();
                _settingsAction.setPauseStatusUpdates(false);
                _settingsAction.getStatusColorAction().setString("C");
            //}
            
        };

    connect(&_settingsAction.getListModelAction(), &VariantAction::variantChanged, this, updateTableModel);

    const auto updateHideShowColumns = [this]() -> void {

        auto shownColumns = _settingsAction.getHiddenShowncolumns().getSelectedOptions();

        QStandardItemModel* model = qobject_cast<QStandardItemModel*>(_settingsAction.getGeneTableView()->model());

        if (model) {
            for (int i = 0; i < model->columnCount(); i++) {
                if (!shownColumns.contains(model->horizontalHeaderItem(i)->text())) {
                    _settingsAction.getGeneTableView()->hideColumn(i);
                }
                else
                {
                    _settingsAction.getGeneTableView()->showColumn(i);

                }
            }
            emit model->layoutChanged();
        }
        };
    connect(&_settingsAction.getHiddenShowncolumns(), &OptionsAction::selectedOptionsChanged, this, updateHideShowColumns);


    const auto triggerInit = [this]() -> void
        {
            int groupIDDeletion = 10;
            int groupId1 = 10 * 2;
            int groupId2 = 10 * 3;
            _settingsAction.removeDatasets(groupIDDeletion);
            _settingsAction.removeDatasets(groupId1);
            _settingsAction.removeDatasets(groupId2);
            _settingsAction.removeDatasets(-1);
            _settingsAction.getMapForHierarchyItemsChangeMethodStopForProjectLoadBlocker().setChecked(false);




            auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
            mv::gui::DatasetPickerAction* colorDatasetPickerAction;
            mv::gui::DatasetPickerAction* pointDatasetPickerAction;
            mv::gui::ViewPluginSamplerAction* samplerActionAction;

            if (scatterplotViewFactory) {
                for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                    if (plugin->getGuiName() == "Scatterplot Embedding View") {
                        pointDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Position"));
                        if (pointDatasetPickerAction) {
                            pointDatasetPickerAction->setCurrentText("");

                            pointDatasetPickerAction->setCurrentDataset(_settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset());

                            colorDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Color"));
                            if (colorDatasetPickerAction)
                            {
                                colorDatasetPickerAction->setCurrentText("");
                                colorDatasetPickerAction->setCurrentDataset(_settingsAction.getBottomClusterNamesDataset().getCurrentDataset());

                            }
                            samplerActionAction = plugin->findChildByPath<mv::gui::ViewPluginSamplerAction>("Sampler");

                            if (samplerActionAction)
                            {
                                samplerActionAction->setTooltipGeneratorFunction([this](const ViewPluginSamplerAction::SampleContext& toolTipContext) -> QString {
                                    QString clusterDatasetId = _settingsAction.getSpeciesNamesDataset().getCurrentDataset().getDatasetId();
                                    return _settingsAction.generateTooltip(toolTipContext, clusterDatasetId, true, "GlobalPointIndices");
                                    });
                            }

                        }
                    }
                }
            }



            

        };

    connect(&mv::projects(), &AbstractProjectManager::projectOpened , this, triggerInit);

    

    //change height of headers

    

    //make long strings in the cells visible and not ...shortened
    _settingsAction.getGeneTableView()->setTextElideMode(Qt::ElideNone);
    _settingsAction.getGeneTableView()->setWordWrap(true);
    _settingsAction.getGeneTableView()->setAlternatingRowColors(false);
    _settingsAction.getGeneTableView()->setSortingEnabled(true);

    //on hovering a cell, show the full text available in a tooltip
    connect(_settingsAction.getGeneTableView(), &QTableView::entered, [this](const QModelIndex& index) {
        if (index.isValid()) {
            QString text = index.data().toString();
            if (!text.isEmpty()) {
                _settingsAction.getGeneTableView()->setToolTip(text);
            }
        }
        });



    _settingsAction.getGeneTableView()->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _settingsAction.getGeneTableView()->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _settingsAction.getGeneTableView()->sortByColumn(1, Qt::DescendingOrder);

    _settingsAction.getGeneTableView()->verticalHeader()->hide();
    _settingsAction.getGeneTableView()->setMouseTracking(true);
    _settingsAction.getGeneTableView()->setToolTipDuration(10000);
    QFont fontLs = _settingsAction.getGeneTableView()->horizontalHeader()->font();
    fontLs.setBold(true);
    _settingsAction.getGeneTableView()->horizontalHeader()->setFont(fontLs);
    _settingsAction.getGeneTableView()->setStyleSheet("QTableView::item:selected { background-color: #00A2ED; }");
    _settingsAction.getGeneTableView()->horizontalHeader()->setHighlightSections(false);
    _settingsAction.getGeneTableView()->verticalHeader()->setHighlightSections(false);






    _settingsAction.getSelectionDetailsTable()->setTextElideMode(Qt::ElideNone);
    _settingsAction.getSelectionDetailsTable()->setWordWrap(true);
    _settingsAction.getSelectionDetailsTable()->setAlternatingRowColors(false);
    _settingsAction.getSelectionDetailsTable()->setSortingEnabled(true);

    //on hovering a cell, show the full text available in a tooltip
    connect(_settingsAction.getSelectionDetailsTable(), &QTableView::entered, [this](const QModelIndex& index) {
        if (index.isValid()) {
            QString text = index.data().toString();
            if (!text.isEmpty()) {
                _settingsAction.getSelectionDetailsTable()->setToolTip(text);
            }
        }
        });



    _settingsAction.getSelectionDetailsTable()->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _settingsAction.getSelectionDetailsTable()->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _settingsAction.getSelectionDetailsTable()->sortByColumn(2, Qt::DescendingOrder);

    _settingsAction.getSelectionDetailsTable()->verticalHeader()->hide();
    _settingsAction.getSelectionDetailsTable()->setMouseTracking(true);
    _settingsAction.getSelectionDetailsTable()->setToolTipDuration(10000);
    QFont fontTbl = _settingsAction.getSelectionDetailsTable()->horizontalHeader()->font();
    fontTbl.setBold(true);
    _settingsAction.getSelectionDetailsTable()->horizontalHeader()->setFont(fontTbl);
    _settingsAction.getSelectionDetailsTable()->setStyleSheet("QTableView::item:selected { background-color: #00A2ED; }");
    _settingsAction.getSelectionDetailsTable()->horizontalHeader()->setHighlightSections(false);
    _settingsAction.getSelectionDetailsTable()->verticalHeader()->setHighlightSections(false);












    auto mainLayout = new QVBoxLayout();
    mainLayout->setContentsMargins(0, 0, 0, 0);
    mainLayout->setSpacing(0);

    auto mainOptionsLayout = new QHBoxLayout();
    mainOptionsLayout->setSpacing(5);
    mainOptionsLayout->setContentsMargins(1, 1, 1, 1);

    auto extraOptionsGroup = new VerticalGroupAction(this, "Settings");
    extraOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("cog"));
    extraOptionsGroup->addAction(&_settingsAction.getListModelAction());
    extraOptionsGroup->addAction(&_settingsAction.getSelectedGeneAction());
    extraOptionsGroup->addAction(&_settingsAction.getSelectedRowIndexAction());
    extraOptionsGroup->addAction(&_settingsAction.getFilteringEditTreeDatasetAction());
    extraOptionsGroup->addAction(&_settingsAction.getOptionSelectionAction());
    extraOptionsGroup->addAction(&_settingsAction.getFilteredGeneNames());
    extraOptionsGroup->addAction(&_settingsAction.getCreateRowMultiSelectTree());
    //extraOptionsGroup->addAction(&_settingsAction.getPerformGeneTableTsneAction());
    extraOptionsGroup->addAction(&_settingsAction.getGeneNamesConnection());
    extraOptionsGroup->addAction(&_settingsAction.getHiddenShowncolumns());
    auto datasetAndLinkerOptionsGroup = new VerticalGroupAction(this, "Dataset and Linker Options");
    datasetAndLinkerOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("link"));
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getReferenceTreeDatasetAction());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getMainPointsDataset());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getEmbeddingDataset());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getSpeciesNamesDataset());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getBottomClusterNamesDataset());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getMiddleClusterNamesDataset());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getTopClusterNamesDataset());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getScatterplotEmbeddingPointsUMAPOption());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getSpeciesExplorerInMap());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getSelctedSpeciesVals());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getStatusColorAction());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getClusterOrderHierarchy());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getRightClickedCluster());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getClearRightClickedCluster());


    auto tsneOptionsGroup = new VerticalGroupAction(this, "Options");
    tsneOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("tools"));
    //tsneOptionsGroup->addAction(&_settingsAction.getUsePreComputedTSNE());
    //tsneOptionsGroup->addAction(&_settingsAction.getTsnePerplexity());
    tsneOptionsGroup->addAction(&_settingsAction.getClusterCountSortingType());
    tsneOptionsGroup->addAction(&_settingsAction.getScatterplotReembedColorOption());
    tsneOptionsGroup->addAction(&_settingsAction.getApplyLogTransformation());
    //tsneOptionsGroup->addAction(&_settingsAction.getTopHierarchyClusterNamesFrequencyInclusionList());
    //tsneOptionsGroup->addAction(&_settingsAction.getMiddleHierarchyClusterNamesFrequencyInclusionList());
    //tsneOptionsGroup->addAction(&_settingsAction.getBottomHierarchyClusterNamesFrequencyInclusionList());
    //tsneOptionsGroup->addAction(&_settingsAction.getComputeTreesToDisplayFromHierarchy());
    //tsneOptionsGroup->addAction(&_settingsAction.getPerformGeneTableTsneAction());
    //tsneOptionsGroup->addAction(&_settingsAction.getPerformGeneTableTsnePerplexity());
    //tsneOptionsGroup->addAction(&_settingsAction.getPerformGeneTableTsneKnn());
    //tsneOptionsGroup->addAction(&_settingsAction.getPerformGeneTableTsneDistance());
    //tsneOptionsGroup->addAction(&_settingsAction.getPerformGeneTableTsneTrigger());
    

    auto mainOptionsGroupLayout = new QVBoxLayout();
    auto mainOptionsGroup1 = new HorizontalGroupAction(this, "MainGroup1");
    auto mainOptionsGroup2 = new HorizontalGroupAction(this, "MainGroup2");
    mainOptionsGroup1->setIcon(Application::getIconFont("FontAwesome").getIcon("database"));
    mainOptionsGroup2->setIcon(Application::getIconFont("FontAwesome").getIcon("play"));

 
    mainOptionsGroup1->addAction(&_settingsAction.getStartComputationTriggerAction());
    mainOptionsGroup1->addAction(&_settingsAction.getTopNGenesFilter());
    mainOptionsGroup1->addAction(&_settingsAction.getTypeofTopNGenes());
    
    
    
    mainOptionsGroup2->addAction(&_settingsAction.getRemoveRowSelection());
    mainOptionsGroup2->addAction(&_settingsAction.getSpeciesExplorerInMapTrigger());
    mainOptionsGroup2->addAction(&_settingsAction.getRevertRowSelectionChangesToInitial());
    mainOptionsGroup2->addAction(&_settingsAction.getToggleScatterplotSelection()); 
    

    auto group1Widget = mainOptionsGroup1->createWidget(&getWidget());
    group1Widget->setMaximumWidth(450);
    mainOptionsGroupLayout->addWidget(group1Widget);

    auto group2Widget = mainOptionsGroup2->createWidget(&getWidget());
    group2Widget->setMaximumWidth(450);
    mainOptionsGroupLayout->addWidget(group2Widget);
    auto statusBarWiddget = _settingsAction.getStatusBarActionWidget();
    statusBarWiddget->setMaximumWidth(600);

    auto searchBoxWidget = _settingsAction.getSearchBox();
    searchBoxWidget->setMaximumWidth(600);

    

    auto startsANDSearchLayout = new QVBoxLayout();
    startsANDSearchLayout->addWidget(statusBarWiddget);
    startsANDSearchLayout->addWidget(searchBoxWidget);
    mainOptionsLayout->addLayout(startsANDSearchLayout);
    mainOptionsLayout->addLayout(mainOptionsGroupLayout);
    //mainOptionsLayout->addWidget(statusBarWiddget);
    //mainOptionsLayout->addLayout(mainOptionsGroupLayout);
    //mainOptionsLayout->addWidget(searchBoxWidget);



    auto linkerandtsneLayout = new QVBoxLayout();
    auto linkerWidget = datasetAndLinkerOptionsGroup->createCollapsedWidget(&getWidget());
    linkerWidget->setMaximumHeight(22);
    //linkerandtsneLayout->addWidget(linkerWidget);
    auto tsneWidget = tsneOptionsGroup->createCollapsedWidget(&getWidget());
    tsneWidget->setMaximumHeight(22);
    linkerandtsneLayout->addWidget(tsneWidget);
    
    mainOptionsLayout->addLayout(linkerandtsneLayout);

    //mainOptionsLayout->addWidget(tsneOptionsGroup->createCollapsedWidget(&getWidget()), 3);
    //mainOptionsLayout->addWidget(datasetAndLinkerOptionsGroup->createCollapsedWidget(&getWidget()), 2);
    //mainOptionsLayout->addWidget(extraOptionsGroup->createCollapsedWidget(&getWidget()), 1);

    auto fullSettingsLayout = new QVBoxLayout();
    fullSettingsLayout->addLayout(mainOptionsLayout);
    mainLayout->addLayout(fullSettingsLayout);

    connect(_settingsAction.getGeneTableView(), &QTableView::entered, [this](const QModelIndex& index) {
        if (index.isValid()) {
            //QString text = index.model()->data(index).toString();
            //QToolTip::showText(QCursor::pos(), text, _settingsAction.getGeneTableView());
            //alternative to above because of linux errors
            QAbstractItemModel* model = _settingsAction.getGeneTableView()->model();
            for (int row = 0; row < model->rowCount(); ++row) {
                for (int column = 0; column < model->columnCount(); ++column) {
                    QModelIndex index = model->index(row, column);
                    QString text = model->data(index).toString();
                    model->setData(index, text, Qt::ToolTipRole);
                }
            }

        }
        });



    QFrame* verticalLine = new QFrame();
    verticalLine->setFrameShape(QFrame::VLine);
    verticalLine->setFrameShadow(QFrame::Sunken);
    verticalLine->setLineWidth(1);

    _settingsAction.getGeneTableView()->setMaximumWidth(165);
    _settingsAction.getTableSplitter()->addWidget(_settingsAction.getGeneTableView());
    _settingsAction.getTableSplitter()->addWidget(verticalLine);
    _settingsAction.getTableSplitter()->addWidget(_settingsAction.getSelectionDetailsTable());

    // Add the splitter to the main layout
    mainLayout->addLayout(_settingsAction.getTableSplitter());

    mainLayout->addLayout(_settingsAction.getSelectedCellClusterInfoStatusBar());
    _settingsAction.getStatusColorAction().setString("M");

    getWidget().setLayout(mainLayout);





}

void CrossSpeciesComparisonGeneDetectPlugin::geneExplorer()
{
    std::vector<std::seed_seq::result_type> selectedSpeciesIndices;
    auto speciesDataset = _settingsAction.getSpeciesNamesDataset().getCurrentDataset();
    auto umapDataset = _settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset();
    auto mainPointsDataset = _settingsAction.getMainPointsDataset().getCurrentDataset();
    std::vector<std::seed_seq::result_type> filtSelectInndx;
    QString gene = _settingsAction.getSelectedGeneAction().getString();
    QStringList finalsettingSpeciesNamesArray = _settingsAction.getSpeciesExplorerInMap().getSelectedOptions();

    if (!speciesDataset.isValid() || !umapDataset.isValid() || !mainPointsDataset.isValid() || !_settingsAction.getFilteredUMAPDatasetPoints().isValid() || !_settingsAction.getFilteredUMAPDatasetColors().isValid() || !_settingsAction.getFilteredUMAPDatasetClusters().isValid())
    {
        qDebug() << "One of the datasets is not valid";
        return;
    }
    if (gene == "")
    {
        qDebug() << "Gene is not selected";
        return;
    }


    if (finalsettingSpeciesNamesArray.isEmpty())
    {
        qDebug() << "Species names are not selected";
        return;
    }


    {
        auto speciesClusterDataset = mv::data().getDataset<Clusters>(speciesDataset.getDatasetId());
        auto umapPointsDataset = mv::data().getDataset<Points>(umapDataset.getDatasetId());
        auto fullMainPointsDataset = mv::data().getDataset<Points>(mainPointsDataset.getDatasetId());

        std::unordered_set<QString> speciesNamesSet(finalsettingSpeciesNamesArray.begin(), finalsettingSpeciesNamesArray.end());
        std::map<QString, std::pair<QColor, std::vector<int>>> selectedFilteredUMAPDatasetColorsMap;
        for (const auto& species : speciesClusterDataset->getClusters()) {
            if (speciesNamesSet.find(species.getName()) != speciesNamesSet.end()) {
                const auto& indices = species.getIndices();
                selectedSpeciesIndices.insert(selectedSpeciesIndices.end(), indices.begin(), indices.end());
                selectedFilteredUMAPDatasetColorsMap[species.getName()] = std::make_pair(species.getColor(), std::vector<int>(indices.begin(), indices.end()));
            }
        }


        std::vector<std::seed_seq::result_type>& selectedIndicesFromStorage = _settingsAction.getSelectedIndicesFromStorage();
        std::unordered_set<std::seed_seq::result_type> indicesSet(selectedIndicesFromStorage.begin(), selectedIndicesFromStorage.end());
        filtSelectInndx.reserve(selectedSpeciesIndices.size());
        for (int i = 0; i < selectedSpeciesIndices.size(); ++i) {
            if (indicesSet.find(selectedSpeciesIndices[i]) != indicesSet.end()) {
                filtSelectInndx.push_back(i);
            }
        }
        auto dimensionNamesUmap = umapPointsDataset->getDimensionNames();
        auto numDimensions = umapPointsDataset->getNumDimensions();
        std::vector<int> geneIndicesSpecies(numDimensions);
        std::iota(geneIndicesSpecies.begin(), geneIndicesSpecies.end(), 0);


        if (selectedSpeciesIndices.size() > 0)
        {
            std::vector<float> resultContainerSpeciesUMAP(selectedSpeciesIndices.size() * umapPointsDataset->getNumDimensions());
            umapPointsDataset->populateDataForDimensions(resultContainerSpeciesUMAP, geneIndicesSpecies, selectedSpeciesIndices);
            auto speciesDataId = _settingsAction.getFilteredUMAPDatasetPoints().getDatasetId();
            int tempnumPoints = static_cast<int>(selectedSpeciesIndices.size());
            int tempNumDimensions = static_cast<int>(geneIndicesSpecies.size());
            _settingsAction.populatePointData(speciesDataId, resultContainerSpeciesUMAP, tempnumPoints, tempNumDimensions, dimensionNamesUmap);



            std::vector<float> resultContainerSpeciesColors(selectedSpeciesIndices.size());
            std::vector<int> selectedGeneIndex;

            auto dimensionNames = fullMainPointsDataset->getDimensionNames();
            auto it = std::find(dimensionNames.begin(), dimensionNames.end(), gene);
            if (it != dimensionNames.end()) {
                selectedGeneIndex.push_back(static_cast<int>(std::distance(dimensionNames.begin(), it)));
            }


            fullMainPointsDataset->populateDataForDimensions(resultContainerSpeciesColors, selectedGeneIndex, selectedSpeciesIndices);
            auto speciesColorDataId = _settingsAction.getFilteredUMAPDatasetColors().getDatasetId();
            auto speciesClusterDataId = _settingsAction.getFilteredUMAPDatasetClusters().getDatasetId();
            int tempnumPointsColors = static_cast<int>(selectedSpeciesIndices.size());

            std::vector<QString> columnGeneColors = { gene };
            int tempNumDimensionsColors = static_cast<int>(columnGeneColors.size());
            if (_settingsAction.getApplyLogTransformation().isChecked())
            {
                applyLogTransformation(resultContainerSpeciesColors);
            }
            _settingsAction.populatePointData(speciesColorDataId, resultContainerSpeciesColors, tempnumPointsColors, tempNumDimensionsColors, columnGeneColors);
            //std::map<QString, std::pair<QColor, std::vector<int>>> selectedFilteredUMAPDatasetColorsMap;
            //TODO: populate the selectedFilteredUMAPDatasetColorsMap

            _settingsAction.populateClusterData(speciesClusterDataId, selectedFilteredUMAPDatasetColorsMap);

            auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
            mv::gui::DatasetPickerAction* colorDatasetPickerAction;
            mv::gui::DatasetPickerAction* pointDatasetPickerAction;
            mv::gui::ViewPluginSamplerAction* samplerActionAction;

            if (scatterplotViewFactory) {
                for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                    if (plugin->getGuiName() == "Scatterplot Embedding View") {
                        //plugin->printChildren();
                        pointDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Position"));
                        if (pointDatasetPickerAction) {
                            pointDatasetPickerAction->setCurrentText("");

                            pointDatasetPickerAction->setCurrentDataset(_settingsAction.getFilteredUMAPDatasetPoints());

                            colorDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Color"));
                            if (colorDatasetPickerAction)
                            {
                                colorDatasetPickerAction->setCurrentText("");
                                colorDatasetPickerAction->setCurrentDataset(_settingsAction.getFilteredUMAPDatasetColors());

                            }

                            auto focusSelectionAction = dynamic_cast<ToggleAction*>(plugin->findChildByPath("Settings/Plot/Point/Focus selection"));

                            if (focusSelectionAction)
                            {
                                focusSelectionAction->setChecked(true);

                            }

                            auto opacityAction = dynamic_cast<DecimalAction*>(plugin->findChildByPath("Settings/Plot/Point/Point opacity/Point opacity"));
                            if (opacityAction)
                            {
                                opacityAction->setValue(20.0);
                            }
                            samplerActionAction = plugin->findChildByPath<mv::gui::ViewPluginSamplerAction>("Sampler");

                            if (samplerActionAction)
                            {
                                samplerActionAction->setTooltipGeneratorFunction([this](const ViewPluginSamplerAction::SampleContext& toolTipContext) -> QString {
                                    QString clusterDatasetId = _settingsAction.getFilteredUMAPDatasetClusters().getDatasetId();
                                    return _settingsAction.generateTooltip(toolTipContext, clusterDatasetId,true, "GlobalPointIndices");
                                    });
                            }
                        }
                    }
                }
            }



            _settingsAction.getFilteredUMAPDatasetPoints()->setSelectionIndices(filtSelectInndx);
            mv::events().notifyDatasetDataSelectionChanged(_settingsAction.getFilteredUMAPDatasetPoints());


        }


    }
}


void CrossSpeciesComparisonGeneDetectPlugin::geneExplorer(QString selectedSpecies)
{
    std::vector<std::seed_seq::result_type> selectedSpeciesIndices;
    auto speciesDataset = _settingsAction.getSpeciesNamesDataset().getCurrentDataset();
    auto umapDataset = _settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset();
    auto mainPointsDataset = _settingsAction.getMainPointsDataset().getCurrentDataset();
    std::vector<std::seed_seq::result_type> filtSelectInndx;
    QString gene = _settingsAction.getSelectedGeneAction().getString();
    //QStringList finalsettingSpeciesNamesArray = _settingsAction.getSpeciesExplorerInMap().getSelectedOptions();

    if (!speciesDataset.isValid() || !umapDataset.isValid() || !mainPointsDataset.isValid() || !_settingsAction.getFilteredUMAPDatasetPoints().isValid() || !_settingsAction.getFilteredUMAPDatasetColors().isValid() || !_settingsAction.getFilteredUMAPDatasetClusters().isValid())
    {
        qDebug() << "One of the datasets is not valid";
        return;
    }
    if (gene == "")
    {
        qDebug() << "Gene is not selected";
        return;
    }


    if (selectedSpecies.isEmpty())
    {
        qDebug() << "Species name empty";
        return;
    }


    {
        auto speciesClusterDataset = mv::data().getDataset<Clusters>(speciesDataset.getDatasetId());
        auto umapPointsDataset = mv::data().getDataset<Points>(umapDataset.getDatasetId());
        auto fullMainPointsDataset = mv::data().getDataset<Points>(mainPointsDataset.getDatasetId());

        std::unordered_set<QString> speciesNamesSet;
        speciesNamesSet.insert(selectedSpecies);
        std::map<QString, std::pair<QColor, std::vector<int>>> selectedFilteredUMAPDatasetColorsMap;
        for (const auto& species : speciesClusterDataset->getClusters()) {
            if (speciesNamesSet.find(species.getName()) != speciesNamesSet.end()) {
                const auto& indices = species.getIndices();
                selectedSpeciesIndices.insert(selectedSpeciesIndices.end(), indices.begin(), indices.end());
                selectedFilteredUMAPDatasetColorsMap[species.getName()] = std::make_pair(species.getColor(), std::vector<int>(indices.begin(), indices.end()));
            }
        }


        std::vector<std::seed_seq::result_type>& selectedIndicesFromStorage = _settingsAction.getSelectedIndicesFromStorage();
        std::unordered_set<std::seed_seq::result_type> indicesSet(selectedIndicesFromStorage.begin(), selectedIndicesFromStorage.end());
        filtSelectInndx.reserve(selectedSpeciesIndices.size());
        for (int i = 0; i < selectedSpeciesIndices.size(); ++i) {
            if (indicesSet.find(selectedSpeciesIndices[i]) != indicesSet.end()) {
                filtSelectInndx.push_back(i);
            }
        }
        auto dimensionNamesUmap = umapPointsDataset->getDimensionNames();
        auto numDimensions = umapPointsDataset->getNumDimensions();
        std::vector<int> geneIndicesSpecies(numDimensions);
        std::iota(geneIndicesSpecies.begin(), geneIndicesSpecies.end(), 0);


        if (selectedSpeciesIndices.size() > 0)
        {
            std::vector<float> resultContainerSpeciesUMAP(selectedSpeciesIndices.size() * umapPointsDataset->getNumDimensions());
            umapPointsDataset->populateDataForDimensions(resultContainerSpeciesUMAP, geneIndicesSpecies, selectedSpeciesIndices);
            auto speciesDataId = _settingsAction.getFilteredUMAPDatasetPoints().getDatasetId();
            int tempnumPoints = static_cast<int>(selectedSpeciesIndices.size());
            int tempNumDimensions = static_cast<int>(geneIndicesSpecies.size());
            _settingsAction.populatePointData(speciesDataId, resultContainerSpeciesUMAP, tempnumPoints, tempNumDimensions, dimensionNamesUmap);



            std::vector<float> resultContainerSpeciesColors(selectedSpeciesIndices.size());
            std::vector<int> selectedGeneIndex;

            auto dimensionNames = fullMainPointsDataset->getDimensionNames();
            auto it = std::find(dimensionNames.begin(), dimensionNames.end(), gene);
            if (it != dimensionNames.end()) {
                selectedGeneIndex.push_back(std::distance(dimensionNames.begin(), it));
            }


            fullMainPointsDataset->populateDataForDimensions(resultContainerSpeciesColors, selectedGeneIndex, selectedSpeciesIndices);
            auto speciesColorDataId = _settingsAction.getFilteredUMAPDatasetColors().getDatasetId();
            auto speciesClusterDataId = _settingsAction.getFilteredUMAPDatasetClusters().getDatasetId();
            int tempnumPointsColors = static_cast<int>(selectedSpeciesIndices.size());

            std::vector<QString> columnGeneColors = { gene };
            int tempNumDimensionsColors = static_cast<int>(columnGeneColors.size());
            if (_settingsAction.getApplyLogTransformation().isChecked())
            {
                applyLogTransformation(resultContainerSpeciesColors);
            }
            _settingsAction.populatePointData(speciesColorDataId, resultContainerSpeciesColors, tempnumPointsColors, tempNumDimensionsColors, columnGeneColors);
            //std::map<QString, std::pair<QColor, std::vector<int>>> selectedFilteredUMAPDatasetColorsMap;
            //TODO: populate the selectedFilteredUMAPDatasetColorsMap

            _settingsAction.populateClusterData(speciesClusterDataId, selectedFilteredUMAPDatasetColorsMap);

            auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
            mv::gui::DatasetPickerAction* colorDatasetPickerAction;
            mv::gui::DatasetPickerAction* pointDatasetPickerAction;
            mv::gui::ViewPluginSamplerAction* samplerActionAction;

            if (scatterplotViewFactory) {
                for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                    if (plugin->getGuiName() == "Scatterplot Embedding View") {
                        //plugin->printChildren();
                        pointDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Position"));
                        if (pointDatasetPickerAction) {
                            pointDatasetPickerAction->setCurrentText("");

                            pointDatasetPickerAction->setCurrentDataset(_settingsAction.getFilteredUMAPDatasetPoints());

                            colorDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Color"));
                            if (colorDatasetPickerAction)
                            {
                                colorDatasetPickerAction->setCurrentText("");
                                colorDatasetPickerAction->setCurrentDataset(_settingsAction.getFilteredUMAPDatasetColors());

                            }

                            auto focusSelectionAction = dynamic_cast<ToggleAction*>(plugin->findChildByPath("Settings/Plot/Point/Focus selection"));

                            if (focusSelectionAction)
                            {
                                focusSelectionAction->setChecked(true);

                            }

                            auto opacityAction = dynamic_cast<DecimalAction*>(plugin->findChildByPath("Settings/Plot/Point/Point opacity/Point opacity"));
                            if (opacityAction)
                            {
                                opacityAction->setValue(20.0);
                            }
                            samplerActionAction = plugin->findChildByPath<mv::gui::ViewPluginSamplerAction>("Sampler");

                            if (samplerActionAction)
                            {
                                samplerActionAction->setTooltipGeneratorFunction([this](const ViewPluginSamplerAction::SampleContext& toolTipContext) -> QString {
                                    QString clusterDatasetId = _settingsAction.getFilteredUMAPDatasetClusters().getDatasetId();
                                    return _settingsAction.generateTooltip(toolTipContext, clusterDatasetId, true, "GlobalPointIndices");
                                    });
                            }
                        }
                    }
                }
            }



            _settingsAction.getFilteredUMAPDatasetPoints()->setSelectionIndices(filtSelectInndx);
            mv::events().notifyDatasetDataSelectionChanged(_settingsAction.getFilteredUMAPDatasetPoints());


        }


    }
}

void CrossSpeciesComparisonGeneDetectPlugin::adjustTableWidths(const QString& value) {
   /*
    // Assuming _settingsAction.getHorizontalLayout() returns your QHBoxLayout
    QHBoxLayout* layout = _settingsAction.getTableSplitter();
    if (!layout) return;

    QWidget* parentWidget = layout->parentWidget();
    if (!parentWidget) return;

    int totalWidth = parentWidget->width();

    double tableViewRatio = 1;
    double selectionDetailsTableRatio = 2;

    if (value == "small") {
        tableViewRatio = 1.2;
        selectionDetailsTableRatio = 1.2;
    }

    int tableViewWidth = static_cast<int>(totalWidth * tableViewRatio / (tableViewRatio + selectionDetailsTableRatio));
    int selectionDetailsTableWidth = totalWidth - tableViewWidth;

    // Assuming the first two widgets in the layout are the ones we want to adjust
    if (layout->count() >= 2) {
        QWidget* tableViewWidget = layout->itemAt(0)->widget();
        QWidget* selectionDetailsTableWidget = layout->itemAt(1)->widget();

        if (tableViewWidget && selectionDetailsTableWidget) {
            tableViewWidget->setMinimumWidth(tableViewWidth);
            tableViewWidget->setMaximumWidth(tableViewWidth);

            selectionDetailsTableWidget->setMinimumWidth(selectionDetailsTableWidth);
            selectionDetailsTableWidget->setMaximumWidth(selectionDetailsTableWidth);
        }
    }
    */
}



QColor getColorFromValue(int value, int min, int max) {
    if (value < min) value = min;
    if (value > max) value = max;

    int range = max - min;
    if (range == 0) return QColor(Qt::gray);

    int blue = 255 * (value - min) / range;

    return QColor(255 - blue, 255 - blue, 255);
}



void CrossSpeciesComparisonGeneDetectPlugin::modifyListData()
{
    try {
    //qDebug() << "It's here";
        auto variant = _settingsAction.getListModelAction().getVariant();
        QStandardItemModel* model = qobject_cast<QStandardItemModel*>(variant.value<QAbstractItemModel*>());

        if (!_settingsAction.getGeneTableView()) {
            // qDebug() << "_settingsAction.getGeneTableView() is null";
            return;
        }

        if (!model) {
            // qDebug() << "Model is null";
            QAbstractItemModel* currentModel = _settingsAction.getGeneTableView()->model();
            if (currentModel) {
                currentModel->removeRows(0, currentModel->rowCount());
                _settingsAction.getGeneTableView()->update();
            }
            else {
                // qDebug() << "TableView model is null";
            }
            return;
        }
        QSortFilterProxyModel* proxyModel = new QSortFilterProxyModel(this);
        proxyModel->setSourceModel(model);
        proxyModel->setFilterCaseSensitivity(Qt::CaseInsensitive);
        proxyModel->setFilterKeyColumn(0);
        _settingsAction.getGeneTableView()->setModel(proxyModel);
        connect(_settingsAction.getSearchBox(), &CustomLineEdit::textboxSelectedForTyping, this, [this, proxyModel, model]() {
            makeAllRowsVisible(_settingsAction.getGeneTableView(),proxyModel);
            });

        connect(_settingsAction.getSearchBox(), &CustomLineEdit::textboxDeselectedNotTypingAnymore, this, [this, proxyModel, model]() {
            proxyModel->setFilterRegularExpression(QRegularExpression());
            updateRowVisibility(_settingsAction.getUniqueReturnGeneList(), _settingsAction.getGeneTableView(),proxyModel);
            });

        connect(_settingsAction.getSearchBox(), &CustomLineEdit::textChanged, this, [this, proxyModel,model](const QString& text) {
            QTimer::singleShot(300, this, [this, proxyModel, text,model]() {
                if (!text.isEmpty()) {
                    //
                    // proxyModel->setFilterWildcard("*" + text + "*");

                    QRegularExpression regExp(text, QRegularExpression::CaseInsensitiveOption);
                    proxyModel->setFilterRegularExpression(regExp);

                }

                });
            });


        _settingsAction.getSearchBox()->setText("");

        updateRowVisibility(_settingsAction.getUniqueReturnGeneList(), _settingsAction.getGeneTableView(),proxyModel);

        // Hide columns not in the shown columns list
        auto shownColumns = _settingsAction.getHiddenShowncolumns().getSelectedOptions();
        for (int i = 0; i < _settingsAction.getGeneTableView()->model()->columnCount(); i++) {
            if (!shownColumns.contains(model->horizontalHeaderItem(i)->text())) {
                _settingsAction.getGeneTableView()->hideColumn(i);
            }
        }

        // Sort and update the table view
        model->sort(1, Qt::DescendingOrder);
        _settingsAction.getGeneTableView()->resizeColumnsToContents();
        _settingsAction.getGeneTableView()->update();


        //disconnect(_settingsAction.getGeneTableView()->selectionModel(), &QItemSelectionModel::currentChanged, this, nullptr);

    connect(_settingsAction.getGeneTableView()->selectionModel(), &QItemSelectionModel::currentChanged, [this](const QModelIndex& current, const QModelIndex& previous) {
        if (!current.isValid()) return;

        QString gene = current.siblingAtColumn(0).data().toString();
        _settingsAction.getSelectedGeneAction().setString(gene);
        _settingsAction.getSelectedRowIndexAction().setString(QString::number(current.row()));
        _settingsAction.getRemoveRowSelection().setEnabled(true);
        //_settingsAction.enableDisableButtonsAutomatically();
        

        std::map<QString, SpeciesDetailsStats> speciesExpressionMap;
        QStringList finalsettingSpeciesNamesArray;
        QString finalSpeciesNameString;
        const auto* model = current.model();
        const int columnCount = model->columnCount();
        const auto initColumnNamesSize = _settingsAction.getInitColumnNames().size(); 
        for (int i = 0; i < columnCount; ++i) {
            const QString columnName = model->headerData(i, Qt::Horizontal).toString();
            const auto data = current.siblingAtColumn(i).data();

            if (columnName == "Gene Appearance Species Names") {
                finalsettingSpeciesNamesArray = data.toString().split(";");
                finalSpeciesNameString = finalsettingSpeciesNamesArray.join(" @%$,$%@ ");
            }
            else if (columnName == "Statistics") {

                if (!finalsettingSpeciesNamesArray.isEmpty()) {
                    std::map<QString, SpeciesDetailsStats> statisticsValues = convertToStatisticsMap(data.toString());
                    speciesExpressionMap = statisticsValues;
                    selectedCellCountStatusBarRemove();
                    selectedCellStatisticsStatusBarAdd(statisticsValues, finalsettingSpeciesNamesArray);
                }
                else {
    
                    qDebug() << "Warning: Gene Appearance Species Names must be processed before Statistics.";
                }
            }
        }



        std::vector<std::seed_seq::result_type> selectedSpeciesIndices;
        auto speciesDataset = _settingsAction.getSpeciesNamesDataset().getCurrentDataset();
        auto umapDataset = _settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset();
        auto mainPointsDataset = _settingsAction.getMainPointsDataset().getCurrentDataset();
        std::vector<std::seed_seq::result_type> filtSelectInndx;


        if (!speciesDataset.isValid() || !umapDataset.isValid() || !mainPointsDataset.isValid() || !_settingsAction.getFilteredUMAPDatasetPoints().isValid() || !_settingsAction.getFilteredUMAPDatasetColors().isValid() || !_settingsAction.getFilteredUMAPDatasetClusters().isValid())
        {
            qDebug() << "One of the datasets is not valid";
            return;
        }


        {
            auto speciesClusterDataset = mv::data().getDataset<Clusters>(speciesDataset.getDatasetId());
            auto umapPointsDataset = mv::data().getDataset<Points>(umapDataset.getDatasetId());
            auto fullMainPointsDataset = mv::data().getDataset<Points>(mainPointsDataset.getDatasetId());
            std::unordered_set<QString> speciesNamesSet(finalsettingSpeciesNamesArray.begin(), finalsettingSpeciesNamesArray.end());
            std::map<QString, std::pair<QColor, std::vector<int>>> selectedFilteredUMAPDatasetColorsMap;
            for (const auto& species : speciesClusterDataset->getClusters()) {
                if (speciesNamesSet.find(species.getName()) != speciesNamesSet.end()) {
                    const auto& indices = species.getIndices();

                    selectedSpeciesIndices.insert(selectedSpeciesIndices.end(), indices.begin(), indices.end());
                    selectedFilteredUMAPDatasetColorsMap[species.getName()] = std::make_pair(species.getColor(), std::vector<int>(indices.begin(), indices.end()));
                }
            }

            
            std::vector<std::seed_seq::result_type>& selectedIndicesFromStorage = _settingsAction.getSelectedIndicesFromStorage();
            std::unordered_set<std::seed_seq::result_type> indicesSet(selectedIndicesFromStorage.begin(), selectedIndicesFromStorage.end());
            filtSelectInndx.reserve(selectedSpeciesIndices.size());
            for (int i = 0; i < selectedSpeciesIndices.size(); ++i) {
                if (indicesSet.find(selectedSpeciesIndices[i]) != indicesSet.end()) {
                    filtSelectInndx.push_back(i);
                }
            }
            auto dimensionNamesUmap = umapPointsDataset->getDimensionNames();
            auto numDimensions = umapPointsDataset->getNumDimensions();
            std::vector<int> geneIndicesSpecies(numDimensions);
            std::iota(geneIndicesSpecies.begin(), geneIndicesSpecies.end(), 0);


            if (selectedSpeciesIndices.size() > 0)
            {
                std::vector<float> resultContainerSpeciesUMAP(selectedSpeciesIndices.size()* umapPointsDataset->getNumDimensions());
                umapPointsDataset->populateDataForDimensions(resultContainerSpeciesUMAP, geneIndicesSpecies, selectedSpeciesIndices);
                auto speciesDataId= _settingsAction.getFilteredUMAPDatasetPoints().getDatasetId();
                int tempnumPoints = static_cast<int>(selectedSpeciesIndices.size());
                int tempNumDimensions = static_cast<int>(geneIndicesSpecies.size());
                _settingsAction.populatePointData(speciesDataId, resultContainerSpeciesUMAP, tempnumPoints, tempNumDimensions, dimensionNamesUmap);



                std::vector<float> resultContainerSpeciesColors(selectedSpeciesIndices.size());
                std::vector<int> selectedGeneIndex;

                auto dimensionNames = fullMainPointsDataset->getDimensionNames();
                auto it = std::find(dimensionNames.begin(), dimensionNames.end(), gene);
                if (it != dimensionNames.end()) {
                    selectedGeneIndex.push_back(std::distance(dimensionNames.begin(), it));
                }


                fullMainPointsDataset->populateDataForDimensions(resultContainerSpeciesColors, selectedGeneIndex, selectedSpeciesIndices);
                auto speciesColorDataId = _settingsAction.getFilteredUMAPDatasetColors().getDatasetId();
                auto speciesClusterDataId = _settingsAction.getFilteredUMAPDatasetClusters().getDatasetId();
                int tempnumPointsColors = static_cast<int>(selectedSpeciesIndices.size());
                
                std::vector<QString> columnGeneColors = { gene };
                int tempNumDimensionsColors = static_cast<int>(columnGeneColors.size());
                if (_settingsAction.getApplyLogTransformation().isChecked())
                {
                    applyLogTransformation(resultContainerSpeciesColors);
                }
                _settingsAction.populatePointData(speciesColorDataId, resultContainerSpeciesColors, tempnumPointsColors, tempNumDimensionsColors, columnGeneColors);

                //std::map<QString, std::pair<QColor, std::vector<int>>> selectedFilteredUMAPDatasetColorsMap;
                //TODO: populate the selectedFilteredUMAPDatasetColorsMap


                _settingsAction.populateClusterData(speciesClusterDataId, selectedFilteredUMAPDatasetColorsMap);


                auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
                mv::gui::DatasetPickerAction* colorDatasetPickerAction;
                mv::gui::DatasetPickerAction* pointDatasetPickerAction;
                mv::gui::ViewPluginSamplerAction* samplerActionAction;

                if (scatterplotViewFactory) {
                    for (auto plugin : mv::plugins().getPluginsByFactory(scatterplotViewFactory)) {
                        if (plugin->getGuiName() == "Scatterplot Embedding View") {
                            //plugin->printChildren();
                            pointDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Position"));
                            if (pointDatasetPickerAction) {
                                pointDatasetPickerAction->setCurrentText("");

                                pointDatasetPickerAction->setCurrentDataset(_settingsAction.getFilteredUMAPDatasetPoints());

                                colorDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Color"));
                                if (colorDatasetPickerAction)
                                {
                                    colorDatasetPickerAction->setCurrentText("");
                                    colorDatasetPickerAction->setCurrentDataset(_settingsAction.getFilteredUMAPDatasetColors());

                                }

                                auto focusSelectionAction = dynamic_cast<ToggleAction*>(plugin->findChildByPath("Settings/Plot/Point/Focus selection"));

                                if (focusSelectionAction)
                                {
                                    focusSelectionAction->setChecked(true);

                                }

                                auto opacityAction = dynamic_cast<DecimalAction*>(plugin->findChildByPath("Settings/Plot/Point/Point opacity/Point opacity"));
                                if (opacityAction)
                                {
                                    opacityAction->setValue(20.0);
                                }
                                
                                samplerActionAction = plugin->findChildByPath<mv::gui::ViewPluginSamplerAction>("Sampler");

                                if (samplerActionAction)
                                {
                                    samplerActionAction->setTooltipGeneratorFunction([this](const ViewPluginSamplerAction::SampleContext& toolTipContext) -> QString {
                                        QString clusterDatasetId = _settingsAction.getFilteredUMAPDatasetClusters().getDatasetId();
                                        return _settingsAction.generateTooltip(toolTipContext, clusterDatasetId,true, "GlobalPointIndices");
                                        });
                                }
                            }
                        }
                    }
                }

                

                _settingsAction.getFilteredUMAPDatasetPoints()->setSelectionIndices(filtSelectInndx);
                mv::events().notifyDatasetDataSelectionChanged(_settingsAction.getFilteredUMAPDatasetPoints());


            }


        }









        auto referenceTreeDataset = _settingsAction.getReferenceTreeDatasetAction().getCurrentDataset();
        if (referenceTreeDataset.isValid()) {
            auto referenceTree = mv::data().getDataset<CrossSpeciesComparisonTree>(referenceTreeDataset.getDatasetId());
            if (referenceTree.isValid()) {
                QJsonObject speciesDataJson = referenceTree->getTreeData();
                updateSpeciesData(speciesDataJson, speciesExpressionMap);
                referenceTree->setTreeData(speciesDataJson);
                events().notifyDatasetDataChanged(referenceTree);
            }
        }

        std::vector<std::seed_seq::result_type> selectedPoints;
        auto speciesColorClusterDataset = _settingsAction.getTsneDatasetSpeciesColors();
        auto tsneDataset = _settingsAction.getSelectedPointsTSNEDataset();
        if (speciesColorClusterDataset.isValid() && tsneDataset.isValid()) {
            for (const auto& species : speciesColorClusterDataset->getClusters()) {
                if (finalsettingSpeciesNamesArray.contains(species.getName())) {
                    const auto& indices = species.getIndices();
                    selectedPoints.insert(selectedPoints.end(), indices.begin(), indices.end());
                }
            }
            //below two lines causing problems TODO: Test
            //tsneDataset->setSelectionIndices(selectedPoints);
            //mv::events().notifyDatasetDataSelectionChanged(tsneDataset);
        }

        if (_settingsAction.getScatterplotReembedColorOption().getCurrentText() == "Expression"&& !_settingsAction.getPerformGeneTableTsneAction().isChecked()) {
            
            
            
            auto expressionColorPointDataset = _settingsAction.getTsneDatasetExpressionColors();
            
            auto selectedPointsMain = _settingsAction.getSelectedPointsDataset();

            if (expressionColorPointDataset.isValid() && selectedPointsMain.isValid()) {

                const int rowSize = expressionColorPointDataset->getNumPoints();

                if (rowSize == selectedPointsMain->getNumPoints())
                {
                    std::vector<float> resultContainerColorPoints(rowSize, -1.0);

                    QString datasetIdEmb = expressionColorPointDataset->getId();

                    std::vector<int> indexOfGene;
                    auto dimsValsTemp = selectedPointsMain->getDimensionNames();
                    auto it = std::find(dimsValsTemp.begin(), dimsValsTemp.end(), gene);
                    if (it != dimsValsTemp.end()) {
                        indexOfGene.push_back(static_cast<int>(it - dimsValsTemp.begin()));
                    }

                    std::vector<int> tempselectIndices(selectedPointsMain->getNumPoints());
                    std::iota(tempselectIndices.begin(), tempselectIndices.end(), 0);

                    if (indexOfGene.size() > 0 && tempselectIndices.size() > 0)
                    {
                        selectedPointsMain->populateDataForDimensions(resultContainerColorPoints, indexOfGene, tempselectIndices);



                        int rowSizeEmbd = rowSize;
                        int columnSizeEmbd = 1;
                        std::vector<QString> columnGeneEmbd = { gene };
                        if (_settingsAction.getApplyLogTransformation().isChecked())
                        {
                            applyLogTransformation(resultContainerColorPoints);
                        }
                        //
                        _settingsAction.populatePointData(datasetIdEmb, resultContainerColorPoints, rowSizeEmbd, columnSizeEmbd, columnGeneEmbd);

                    }
                }



            }
        }
        else
        {
            //TODO: Add selected gene point selection to scatterplot
            unsigned int selectedGeneIndex = -1;
            if (_settingsAction.getGeneOrder().size() > 0 && _settingsAction.getSelectedPointsTSNEDatasetForGeneTable().isValid())
            {
                auto it = std::find(_settingsAction.getGeneOrder().begin(), _settingsAction.getGeneOrder().end(), gene);
                if (it != _settingsAction.getGeneOrder().end()) {
                    selectedGeneIndex = it - _settingsAction.getGeneOrder().begin();
                }
                //gene
                unsigned int maxIndex = _settingsAction.getSelectedPointsTSNEDatasetForGeneTable()->getNumPoints();
                if (selectedGeneIndex <= maxIndex)
                {
                    _settingsAction.getSelectedPointsTSNEDatasetForGeneTable()->setSelectionIndices({ selectedGeneIndex });
                    mv::events().notifyDatasetDataSelectionChanged(_settingsAction.getSelectedPointsTSNEDatasetForGeneTable());

                }

            }
            


        }

        if (_settingsAction.getSelctedSpeciesVals().getString() == finalSpeciesNameString)
        {
            _settingsAction.getSelctedSpeciesVals().setString("");
        }
        _settingsAction.getSelctedSpeciesVals().setString(finalSpeciesNameString);

        });

    emit model->layoutChanged();

    

    selectedCellStatisticsStatusBarRemove();
    selectedCellCountStatusBarAdd();

        }
        catch (const std::exception& e) {
            qDebug() << "An exception occurred in modifyListData: " << e.what();
        }
        catch (...) {
            qDebug() << "An unknown exception occurred in modifyListData";
        }



}
void CrossSpeciesComparisonGeneDetectPlugin::selectedCellCountStatusBarAdd()
{
    if (!_settingsAction.getSelectedSpeciesCellCountMap().empty())
    {
        // Create a new model for the table view
        QStandardItemModel* model = new QStandardItemModel();

        QSet<QString> topSet = _settingsAction.getCurrentHierarchyItemsTopForTable();
        QSet<QString> middleSet = _settingsAction.getCurrentHierarchyItemsMiddleForTable();
        bool singleColumn = false;
        bool isSameStringValues = false; // Ensure this variable is declared

        //qDebug() << "Top Set:" << topSet;
        //qDebug() << "Middle Set:" << middleSet;

        if (topSet == middleSet) {
            isSameStringValues = true;
        } else {
            isSameStringValues = false;
        }

        QString headerStringToAdd = "";
        if (!isSameStringValues)
        {
            if (middleSet.size() > 1)
            {
                headerStringToAdd = "Combined";
                singleColumn = false;
            }
            else if (middleSet.size() == 1)
            {
                headerStringToAdd = *middleSet.begin();
                singleColumn = false;
            }
        }
        else
        {
            headerStringToAdd = "Abundance";
            singleColumn = true;
        }

        model->setHorizontalHeaderLabels({ "Species", "Fraction of Neuronal", "Fraction of " + headerStringToAdd , "Count Selected", "Count All" });

        for (const auto& [species, details] : _settingsAction.getSelectedSpeciesCellCountMap()) {
            QColor backgroundColor = QColor(details.color); // Ensure color is converted to QColor
            qreal brightness = backgroundColor.lightnessF();
            QColor textColor = (brightness > 0.4) ? Qt::black : Qt::white;

            QList<QStandardItem*> rowItems;

            // Species column
            QStandardItem* speciesItem = new QStandardItem(species);
            speciesItem->setBackground(backgroundColor);
            speciesItem->setForeground(textColor);
            rowItems << speciesItem;

            // Fraction of Neuronal column
            QStandardItem* item = new QStandardItem();
            float topAbundance = 0.0;
            if (details.abundanceTop != 0)
            {
                topAbundance = (static_cast<float>(details.countAbundanceNumerator) / static_cast<float>(details.abundanceTop)) * 100;
            }
            QString formattedValueTop = QString::number(topAbundance, 'f', 2);
            item->setData(QVariant(formattedValueTop), Qt::EditRole);
            rowItems << item;

            // Fraction of Middle column
            item = new QStandardItem();
            float middleAbundance = 0.0;
            if (details.abundanceMiddle != 0)
            {
                middleAbundance = (static_cast<float>(details.countAbundanceNumerator) / static_cast<float>(details.abundanceMiddle)) * 100;
            }
            QString formattedValueMiddle = QString::number(middleAbundance, 'f', 2);
            item->setData(QVariant(formattedValueMiddle), Qt::EditRole);
            rowItems << item;

            // Count Selected column
            item = new QStandardItem();
            item->setData(QVariant(details.selectedCellsCount), Qt::EditRole);
            rowItems << item;

            // Count All column
            item = new QStandardItem();
            auto total = details.selectedCellsCount + details.nonSelectedCellsCount;
            item->setData(QVariant(total), Qt::EditRole);
            rowItems << item;

            model->appendRow(rowItems);
        }

        // Set the model to the table view
        _settingsAction.getSelectionDetailsTable()->setModel(model);

        // Sort the model by the fourth column (Count Selected) in descending order
        model->sort(3, Qt::DescendingOrder);

        // Configure the table view
        _settingsAction.getSelectionDetailsTable()->setSelectionMode(QAbstractItemView::NoSelection);
        _settingsAction.getSelectionDetailsTable()->verticalHeader()->hide();
        if (singleColumn)
        {
            _settingsAction.getSelectionDetailsTable()->hideColumn(2);
        }
        _settingsAction.getSelectionDetailsTable()->resizeColumnsToContents();
        _settingsAction.getSelectionDetailsTable()->update();

        // Emit layoutChanged signal to notify views of the model change
        emit model->layoutChanged();
    }

    adjustTableWidths("small");
}



void CrossSpeciesComparisonGeneDetectPlugin::selectedCellStatisticsStatusBarAdd(std::map<QString, SpeciesDetailsStats> statisticsValues, QStringList selectedSpecies)
{
    if (!_settingsAction.getSelectedSpeciesCellCountMap().empty())
    {
        // Create a new model for the table view
        QStandardItemModel* model = new QStandardItemModel();

        QSet<QString> topSet = _settingsAction.getCurrentHierarchyItemsTopForTable();
        QSet<QString> middleSet = _settingsAction.getCurrentHierarchyItemsMiddleForTable();
        bool singleColumn = false;
        bool isSameStringValues = false; // Ensure this variable is declared

        //qDebug() << "Top Set:" << topSet;
        //qDebug() << "Middle Set:" << middleSet;

        if (topSet == middleSet) {
            isSameStringValues = true;
        }
        else {
            isSameStringValues = false;
        }

        QString headerStringToAdd = "";
        if (!isSameStringValues)
        {
            if (middleSet.size() > 1)
            {
                headerStringToAdd = "Combined";
                singleColumn = false;
            }
            else if (middleSet.size() == 1)
            {
                headerStringToAdd = *middleSet.begin();
                singleColumn = false;
            }
        }
        else
        {
            headerStringToAdd = "Abundance";
            singleColumn = true;
        }

        model->setHorizontalHeaderLabels({ "Species", "Mean Difference", "Appearance Rank", "Fraction of Neuronal", "Fraction of " + headerStringToAdd, "Count Selected", "Mean Selected", "Count Non Selected", "Mean Non Selected" });

        auto colorValues = _settingsAction.getSystemModeColor();
        auto systemColor = colorValues[0];
        auto valuesColor = colorValues[1];

        // Populate the model with sorted data and statistics
        for (const auto& [species, details] : _settingsAction.getSelectedSpeciesCellCountMap()) {
            QColor backgroundColor = QColor(details.color);

            // Calculate the brightness of the background color
            qreal brightness = backgroundColor.lightnessF();

            // Choose text color based on the brightness of the background color
            QColor textColor = (brightness > 0.4) ? Qt::black : Qt::white;

            QList<QStandardItem*> rowItems;
            QString speciesCopy = species;
            speciesCopy.replace("_", " ");
            QStandardItem* item = new QStandardItem(species);
            item->setBackground(backgroundColor);
            item->setForeground(textColor);
            rowItems << item;

            // Find statistics for the species
            auto it = statisticsValues.find(species);
            if (it != statisticsValues.end()) {
                QStandardItem* item;

                item = new QStandardItem();
                float difference = (it->second.meanSelected - it->second.meanNonSelected);
                item->setData(QVariant(QString::number(difference, 'f', 2)), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.rank), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                float topAbundance = 0.0;
                if (it->second.abundanceTop != 0) {
                    topAbundance = (static_cast<float>(it->second.countAbundanceNumerator) / static_cast<float>(it->second.abundanceTop)) * 100;
                }
                QString formattedValueTop = QString::number(topAbundance, 'f', 2);
                item->setData(QVariant(formattedValueTop), Qt::EditRole);
                rowItems << item;

                float middleAbundance = 0.0;
                if (it->second.abundanceMiddle != 0)
                {
                    middleAbundance = (static_cast<float>(it->second.countAbundanceNumerator) / static_cast<float>(it->second.abundanceMiddle)) * 100;
                }
                item = new QStandardItem();
                QString formattedValueMiddle = QString::number(middleAbundance, 'f', 2);
                item->setData(QVariant(formattedValueMiddle), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.countSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                QString formattedMeanSelected = QString::number(it->second.meanSelected, 'f', 2);
                item->setData(QVariant(formattedMeanSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.countNonSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                QString formattedMeanNonSelected = QString::number(it->second.meanNonSelected, 'f', 2);
                item->setData(QVariant(formattedMeanNonSelected), Qt::EditRole);
                rowItems << item;
            }
            else {
                qDebug() << "Species: " + species + " not found in map auto it = statisticsValues.find(species);";
                for (auto it = statisticsValues.begin(); it != statisticsValues.end(); ++it)
                {
                    qDebug() << it->first;
                }
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
            }

            if (selectedSpecies.contains(species)) {
                for (int i = 1; i < rowItems.size(); ++i) {
                    rowItems[i]->setBackground(QBrush(QColor("#00A2ED")));
                    rowItems[i]->setForeground(QBrush(QColor(valuesColor)));
                }
            }
            else
            {
                for (int i = 1; i < rowItems.size(); ++i) {
                    rowItems[i]->setBackground(QBrush(QColor(systemColor)));
                    rowItems[i]->setForeground(QBrush(QColor(valuesColor)));
                }
            }
            _settingsAction.getSpeciesExplorerInMap().setSelectedOptions(selectedSpecies);

            model->appendRow(rowItems);
        }
        _settingsAction.getSelectionDetailsTable()->setModel(model);
        model->sort(2, _settingsAction.getTypeofTopNGenes().getCurrentText() == "Positive" || _settingsAction.getTypeofTopNGenes().getCurrentText() == "Absolute" ? Qt::AscendingOrder : Qt::DescendingOrder);

        _settingsAction.getSelectionDetailsTable()->setSelectionMode(QAbstractItemView::SingleSelection);
        _settingsAction.getSelectionDetailsTable()->verticalHeader()->hide();
        if (singleColumn) {
            _settingsAction.getSelectionDetailsTable()->hideColumn(4);
        }
        _settingsAction.getSelectionDetailsTable()->resizeColumnsToContents();
        _settingsAction.getSelectionDetailsTable()->update();



        emit model->layoutChanged();

        // Add a connect for row selection
        connect(_settingsAction.getSelectionDetailsTable()->selectionModel(), &QItemSelectionModel::selectionChanged, [this](const QItemSelection& selected, const QItemSelection& deselected) {
            static QModelIndex lastSelectedIndex;

            if (selected.isEmpty()) {
                return;
            }
            QModelIndexList selectedIndexes = selected.indexes();
            if (selectedIndexes.isEmpty()) {
                return;
            }
            QModelIndex selectedIndex = selectedIndexes.at(0);
            if (!selectedIndex.isValid()) {
                return;
            }

            lastSelectedIndex = selectedIndex;
            QString species = selectedIndex.siblingAtColumn(0).data().toString();
            QStringList autoSpecies = { species };

            if (!(autoSpecies.size() == 1 && autoSpecies.first().isEmpty())) {
                if (_settingsAction.getSpeciesExplorerInMap().getNumberOfOptions() > 0)
                {
                    _settingsAction.getSpeciesExplorerInMap().setSelectedOptions(autoSpecies);
                    geneExplorer();
                }
            }
            _settingsAction.getClearRightClickedCluster().trigger();
            });
    }
    adjustTableWidths("large");
}



void CrossSpeciesComparisonGeneDetectPlugin::selectedCellCountStatusBarRemove()
{
    _settingsAction.getSelectionDetailsTable()->setModel(new QStandardItemModel());
}

void CrossSpeciesComparisonGeneDetectPlugin::selectedCellStatisticsStatusBarRemove()
{
    _settingsAction.getSelectionDetailsTable()->setModel(new QStandardItemModel());
}


void CrossSpeciesComparisonGeneDetectPlugin::updateSpeciesData(QJsonObject& node, const std::map<QString, SpeciesDetailsStats>& speciesExpressionMap) {
    // Check if the "name" key exists in the current node
    if (node.contains("name")) {
        QString nodeName = node["name"].toString();
        auto it = speciesExpressionMap.find(nodeName);

        if (it != speciesExpressionMap.end()) {
            node["mean"] = std::round(it->second.meanSelected * 100.0) / 100.0; // Round to 2 decimal places
            float differentialVal = it->second.meanSelected - it->second.meanNonSelected;
            node["differential"] = std::round(differentialVal * 100.0) / 100.0; // Round to 2 decimal places
            int rank = it->second.rank;
            node["rank"] = rank;
            node["gene"] = _settingsAction.getSelectedGeneAction().getString();



            float topAbundance = 0.0;
            if (it->second.abundanceTop != 0)
            {
                topAbundance =  (static_cast<float>(it->second.countAbundanceNumerator) / static_cast<float>(it->second.abundanceTop)) * 100;
            }

            float middleAbundance = 0.0;
            if (it->second.abundanceMiddle != 0)
            {
                middleAbundance = (static_cast<float>(it->second.countAbundanceNumerator) / static_cast<float>(it->second.abundanceMiddle)) * 100;
            }

            node["abundanceTop"] = std::round(topAbundance * 100.0) / 100.0;
            node["abundanceMiddle"] =  std::round(middleAbundance * 100.0) / 100.0;
        }
        if (it != speciesExpressionMap.end()) {
            node["cellCounts"] = it->second.countSelected;

        }
    }

    // If the node has "children", recursively update them as well
    if (node.contains("children")) {
        QJsonArray children = node["children"].toArray();
        for (int i = 0; i < children.size(); ++i) {
            QJsonObject child = children[i].toObject();
            updateSpeciesData(child, speciesExpressionMap); // Recursive call
            children[i] = child; // Update the modified object back into the array
        }
        node["children"] = children; // Update the modified array back into the parent JSON object
    }
}
/*
void CrossSpeciesComparisonGeneDetectPlugin::onDataEvent(mv::DatasetEvent* dataEvent)
{
    // Get smart pointer to dataset that changed
    const auto changedDataSet = dataEvent->getDataset();

    // Get GUI name of the dataset that changed
    const auto datasetGuiName = changedDataSet->getGuiName();

    // The data event has a type so that we know what type of data event occurred (e.g. data added, changed, removed, renamed, selection changes)
    switch (dataEvent->getType()) {

        // A points dataset was added
        case EventType::DatasetAdded:
        {
            // Cast the data event to a data added event
            const auto dataAddedEvent = static_cast<DatasetAddedEvent*>(dataEvent);

            // Get the GUI name of the added points dataset and print to the console
            qDebug() << datasetGuiName << "was added";

            break;
        }

        // Points dataset data has changed
        case EventType::DatasetDataChanged:
        {
            // Cast the data event to a data changed event
            const auto dataChangedEvent = static_cast<DatasetDataChangedEvent*>(dataEvent);

            // Get the name of the points dataset of which the data changed and print to the console
            qDebug() << datasetGuiName << "data changed";

            break;
        }

        // Points dataset data was removed
        case EventType::DatasetRemoved:
        {
            // Cast the data event to a data removed event
            const auto dataRemovedEvent = static_cast<DatasetRemovedEvent*>(dataEvent);

            // Get the name of the removed points dataset and print to the console
            qDebug() << datasetGuiName << "was removed";

            break;
        }

        // Points dataset selection has changed
        case EventType::DatasetDataSelectionChanged:
        {
            // Cast the data event to a data selection changed event
            const auto dataSelectionChangedEvent = static_cast<DatasetDataSelectionChangedEvent*>(dataEvent);

            // Get the selection set that changed
            const auto& selectionSet = changedDataSet->getSelection<Points>();

            // Print to the console
            qDebug() << datasetGuiName << "selection has changed";

            break;
        }

        default:
            break;
    }
}
*/

void CrossSpeciesComparisonGeneDetectPlugin::fromVariantMap(const QVariantMap& variantMap)
{
    ViewPlugin::fromVariantMap(variantMap);

    mv::util::variantMapMustContain(variantMap, "CSCGDV:CrossSpeciesComparison Gene Detect Plugin Settings");
    _settingsAction.fromVariantMap(variantMap["CSCGDV:CrossSpeciesComparison Gene Detect Plugin Settings"].toMap());
   // modifyTableData();
   // _settingsAction.getStartComputationTriggerAction().trigger();

}

QVariantMap CrossSpeciesComparisonGeneDetectPlugin::toVariantMap() const
{
    QVariantMap variantMap = ViewPlugin::toVariantMap();

    _settingsAction.insertIntoVariantMap(variantMap);

    return variantMap;
}
ViewPlugin* CrossSpeciesComparisonGeneDetectPluginFactory::produce()
{
    return new CrossSpeciesComparisonGeneDetectPlugin(this);
}

mv::DataTypes CrossSpeciesComparisonGeneDetectPluginFactory::supportedDataTypes() const
{
    DataTypes supportedTypes;

    // This example analysis plugin is compatible with points datasets
    supportedTypes.append(PointType);

    return supportedTypes;
}

mv::gui::PluginTriggerActions CrossSpeciesComparisonGeneDetectPluginFactory::getPluginTriggerActions(const mv::Datasets& datasets) const
{
    PluginTriggerActions pluginTriggerActions;
    /*
    const auto getPluginInstance = [this]() -> CrossSpeciesComparisonGeneDetectPlugin* {
        return dynamic_cast<CrossSpeciesComparisonGeneDetectPlugin*>(plugins().requestViewPlugin(getKind()));
    };

    const auto numberOfDatasets = datasets.count();

    if (numberOfDatasets >= 1 && PluginFactory::areAllDatasetsOfTheSameType(datasets, PointType)) {
        auto pluginTriggerAction = new PluginTriggerAction(const_cast<CrossSpeciesComparisonGeneDetectPluginFactory*>(this), this, "CrossSpeciesComparisonGeneDetect View", "View gene data", getIcon(), [this, getPluginInstance, datasets](PluginTriggerAction& pluginTriggerAction) -> void {
            for (auto dataset : datasets)
                getPluginInstance();
        });

        pluginTriggerActions << pluginTriggerAction;
    }
    */
    return pluginTriggerActions;
}

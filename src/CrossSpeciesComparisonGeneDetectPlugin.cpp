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
#include <execution>
#include<QTooltip>
Q_PLUGIN_METADATA(IID "studio.manivault.CrossSpeciesComparisonGeneDetectPlugin")

using namespace mv;


void applyLogTransformation(std::vector<float>& values) {
    std::transform(std::execution::par, values.begin(), values.end(), values.begin(),
        [](float value) { return std::log(value + 1); });
}

std::map<QString, Statistics> convertToStatisticsMap(const QString& formattedStatistics) {
    std::map<QString, Statistics> statisticsMap;

    // Split the QString into individual species statistics using ";" as a delimiter
    // For compatibility with Qt versions before 5.14, use QString::SplitBehavior enum
    QStringList speciesStatsList = formattedStatistics.split(";", Qt::SkipEmptyParts); // Qt 5.14 and later
    // For Qt versions before 5.14, use the following line instead:
    // QStringList speciesStatsList = formattedStatistics.split(";", QString::SkipEmptyParts); // Before Qt 5.14

    // Regular expression to match the pattern of each statistic
    QRegularExpression regex("Species: (.*), MeanSelected: ([\\d.]+), MedianSelected: ([\\d.]+), ModeSelected: ([\\d.]+), RangeSelected: ([\\d.]+), CountSelected: (\\d+), MeanNotSelected: ([\\d.]+), MedianNotSelected: ([\\d.]+), ModeNotSelected: ([\\d.]+), RangeNotSelected: ([\\d.]+), CountNotSelected: (\\d+)");


    for (const QString& speciesStats : speciesStatsList) {
        QRegularExpressionMatch match = regex.match(speciesStats.trimmed());
        if (match.hasMatch()) {
            QString species = match.captured(1);
            Statistics stats = {
                match.captured(2).toFloat(),
                match.captured(3).toFloat(),
                match.captured(4).toFloat(),
                match.captured(5).toFloat(),
                match.captured(6).toInt(),
                match.captured(7).toFloat(),
                match.captured(8).toFloat(),
                match.captured(9).toFloat(),
                match.captured(10).toFloat(),
                match.captured(11).toInt()
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
                    if (treeDataset.isValid() && _settingsAction.getTableView() && selectedRow >= 0)
                    {
                        QString treeData = _settingsAction.getTableView()->model()->index(selectedRow, 2).data().toString();
                        //qDebug()<< "Tree data: " << treeData;
                        if (!treeData.isEmpty())
                        {

                            QJsonObject valueStringReference = QJsonDocument::fromJson(treeData.toUtf8()).object();
                            if (!valueStringReference.isEmpty())
                            {
                                treeDataset->setTreeData(valueStringReference);
                                events().notifyDatasetDataChanged(treeDataset);
                                //QString firstColumnValue = _settingsAction.getTableView()->model()->index(selectedRow, 0).data().toString();
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
                    firstColumnValues << _settingsAction.getTableView()->model()->index(row, 0).data().toString();
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
            if (_settingsAction.getTableView() && _settingsAction.getTableView()->selectionModel()) {
                // Clear the current index if there's no selection
                _settingsAction.getTableView()->clearSelection();

                // Temporarily disable the selection mode to remove highlight
                QAbstractItemView::SelectionMode oldMode = _settingsAction.getTableView()->selectionMode();
                _settingsAction.getTableView()->setSelectionMode(QAbstractItemView::NoSelection);

                // Clear the current index
                _settingsAction.getTableView()->selectionModel()->setCurrentIndex(QModelIndex(), QItemSelectionModel::NoUpdate);

                // Restore the original selection mode
                _settingsAction.getTableView()->setSelectionMode(oldMode);
                // Update the view to ensure changes are reflected
                _settingsAction.getTableView()->update();
                _settingsAction.getSelctedSpeciesVals().setString("");


                if (_settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset().isValid() && _settingsAction.getClusterNamesDataset().getCurrentDataset().isValid())
                {

                    auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
                    mv::gui::DatasetPickerAction* colorDatasetPickerAction;
                    mv::gui::DatasetPickerAction* pointDatasetPickerAction;


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
                                        colorDatasetPickerAction->setCurrentDataset(_settingsAction.getClusterNamesDataset().getCurrentDataset());

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

            _settingsAction.getRemoveRowSelection().setDisabled(true);
            _settingsAction.getStatusColorAction().setString(statusString);
            selectedCellStatisticsStatusBarRemove();
            selectedCellCountStatusBarAdd();

        };

    connect(&_settingsAction.getRemoveRowSelection(), &TriggerAction::triggered, this, removeRowSelectionTable);

    const auto updateTableModel = [this]() -> void
        {
            modifyTableData();
            _settingsAction.getStatusColorAction().setString("C");
        };

    connect(&_settingsAction.getTableModelAction(), &VariantAction::variantChanged, this, updateTableModel);

    const auto updateHideShowColumns = [this]() -> void {

        auto shownColumns = _settingsAction.getHiddenShowncolumns().getSelectedOptions();

        QStandardItemModel* model = qobject_cast<QStandardItemModel*>(_settingsAction.getTableView()->model());

        if (model) {
            for (int i = 0; i < model->columnCount(); i++) {
                if (!shownColumns.contains(model->horizontalHeaderItem(i)->text())) {
                    _settingsAction.getTableView()->hideColumn(i);
                }
                else
                {
                    _settingsAction.getTableView()->showColumn(i);

                }
            }
            emit model->layoutChanged();
        }
        };
    connect(&_settingsAction.getHiddenShowncolumns(), &OptionsAction::selectedOptionsChanged, this, updateHideShowColumns);



    //change height of headers

    

    //make long strings in the cells visible and not ...shortened
    //_settingsAction.getTableView()->setTextElideMode(Qt::ElideNone);
    //_settingsAction.getTableView()->setWordWrap(true);
    //_settingsAction.getTableView()->setAlternatingRowColors(true);
    //_settingsAction.getTableView()->setSortingEnabled(true);

    //on hovering a cell, show the full text available in a tooltip
    connect(_settingsAction.getTableView(), &QTableView::entered, [this](const QModelIndex& index) {
        if (index.isValid()) {
            QString text = index.data().toString();
            if (!text.isEmpty()) {
                _settingsAction.getTableView()->setToolTip(text);
            }
        }
        });


    /*
    connect(_settingsAction.getTableView(), &QTableView::clicked, [this](const QModelIndex& index) {
        QModelIndex firstColumnIndex = index.sibling(index.row(), 0);
        auto gene = firstColumnIndex.data().toString();
        _settingsAction.getSelectedGeneAction().setString(gene);

        //if (QApplication::keyboardModifiers() & Qt::ShiftModifier) 
        
        //{
            // If Shift is pressed, add the row to the selection
          //  _settingsAction.getTableView()->selectionModel()->select(index, QItemSelectionModel::Select | QItemSelectionModel::Rows);
        //}
        //else {
            // If Shift is not pressed, select only this row
            //_settingsAction.getTableView()->selectionModel()->clearSelection();
           // _settingsAction.getTableView()->selectionModel()->select(index, QItemSelectionModel::Select | QItemSelectionModel::Rows);
       // }

        // Get the selected rows and convert them to a string list
        //QModelIndexList selectedRows = _settingsAction.getTableView()->selectionModel()->selectedRows();
        QStringList selectedRowsStrList;
        for (const QModelIndex& selectedIndex : selectedRows) {
            selectedRowsStrList << QString::number(selectedIndex.row());
        }

        // Join the string list into a single string with comma separation
        QString selectedRowsStr = selectedRowsStrList.join(",");
        _settingsAction.getSelectedRowIndexAction().setString(selectedRowsStr);
        });
   */


    _settingsAction.getTableView()->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _settingsAction.getTableView()->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _settingsAction.getTableView()->sortByColumn(3, Qt::DescendingOrder);

    _settingsAction.getTableView()->verticalHeader()->hide();
    _settingsAction.getTableView()->setMouseTracking(true);
    _settingsAction.getTableView()->setToolTipDuration(10000);
    QFont font = _settingsAction.getTableView()->horizontalHeader()->font();
    font.setBold(true);
    _settingsAction.getTableView()->horizontalHeader()->setFont(font);
    _settingsAction.getTableView()->setStyleSheet("QTableView::item:selected { background-color: #00A2ED; }");
    _settingsAction.getTableView()->horizontalHeader()->setHighlightSections(false);
    _settingsAction.getTableView()->verticalHeader()->setHighlightSections(false);


    auto mainLayout = new QVBoxLayout();
    mainLayout->setContentsMargins(0, 0, 0, 0);
    mainLayout->setSpacing(0);

    auto mainOptionsLayout = new QHBoxLayout();
    mainOptionsLayout->setSpacing(0);
    mainOptionsLayout->setContentsMargins(0, 0, 0, 0);

    auto extraOptionsGroup = new VerticalGroupAction(this, "Settings");
    extraOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("cog"));
    extraOptionsGroup->addAction(&_settingsAction.getTableModelAction());
    extraOptionsGroup->addAction(&_settingsAction.getSelectedGeneAction());
    extraOptionsGroup->addAction(&_settingsAction.getSelectedRowIndexAction());
    extraOptionsGroup->addAction(&_settingsAction.getFilteringEditTreeDatasetAction());
    extraOptionsGroup->addAction(&_settingsAction.getOptionSelectionAction());
    extraOptionsGroup->addAction(&_settingsAction.getFilteredGeneNames());
    extraOptionsGroup->addAction(&_settingsAction.getCreateRowMultiSelectTree());
    extraOptionsGroup->addAction(&_settingsAction.getPerformGeneTableTsneAction());

    auto datasetAndLinkerOptionsGroup = new VerticalGroupAction(this, "Dataset and Linker Options");
    datasetAndLinkerOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("link"));
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getReferenceTreeDatasetAction());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getMainPointsDataset());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getEmbeddingDataset());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getSpeciesNamesDataset());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getClusterNamesDataset());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getScatterplotEmbeddingPointsUMAPOption());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getGeneNamesConnection());
    datasetAndLinkerOptionsGroup->addAction(&_settingsAction.getSelctedSpeciesVals());

    auto tsneOptionsGroup = new VerticalGroupAction(this, "Options");
    tsneOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("tools"));
    tsneOptionsGroup->addAction(&_settingsAction.getUsePreComputedTSNE());
    tsneOptionsGroup->addAction(&_settingsAction.getTsnePerplexity());
    tsneOptionsGroup->addAction(&_settingsAction.getTypeofTopNGenes());
    tsneOptionsGroup->addAction(&_settingsAction.getHiddenShowncolumns());

    auto mainOptionsGroupLayout = new QVBoxLayout();
    auto mainOptionsGroup1 = new HorizontalGroupAction(this, "MainGroup1");
    auto mainOptionsGroup2 = new HorizontalGroupAction(this, "MainGroup2");
    mainOptionsGroup1->setIcon(Application::getIconFont("FontAwesome").getIcon("database"));
    mainOptionsGroup2->setIcon(Application::getIconFont("FontAwesome").getIcon("play"));
    mainOptionsGroup1->addAction(&_settingsAction.getTopNGenesFilter());
    mainOptionsGroup1->addAction(&_settingsAction.getScatterplotReembedColorOption());

    mainOptionsGroup2->addAction(&_settingsAction.getStartComputationTriggerAction());
    mainOptionsGroup2->addAction(&_settingsAction.getRemoveRowSelection());
    

    auto group1Widget = mainOptionsGroup1->createWidget(&getWidget());
    group1Widget->setMaximumWidth(460);
    mainOptionsGroupLayout->addWidget(group1Widget);

    auto group2Widget = mainOptionsGroup2->createWidget(&getWidget());
    group2Widget->setMaximumWidth(500);
    mainOptionsGroupLayout->addWidget(group2Widget);

    mainOptionsLayout->addWidget(_settingsAction.getStatusBarActionWidget());
    mainOptionsLayout->addLayout(mainOptionsGroupLayout);
    mainOptionsLayout->addWidget(tsneOptionsGroup->createCollapsedWidget(&getWidget()), 3);
    mainOptionsLayout->addWidget(datasetAndLinkerOptionsGroup->createCollapsedWidget(&getWidget()), 2);
    mainOptionsLayout->addWidget(extraOptionsGroup->createCollapsedWidget(&getWidget()), 1);

    auto fullSettingsLayout = new QVBoxLayout();
    fullSettingsLayout->addLayout(mainOptionsLayout);
    mainLayout->addLayout(fullSettingsLayout);

    connect(_settingsAction.getTableView(), &QTableView::entered, [this](const QModelIndex& index) {
        if (index.isValid()) {
            QString text = index.model()->data(index).toString();
            QToolTip::showText(QCursor::pos(), text, _settingsAction.getTableView());
        }
        });

    _settingsAction.getTableSplitter()->addWidget(_settingsAction.getTableView(),1.2);
    _settingsAction.getTableSplitter()->addWidget(_settingsAction.getSelectionDetailsTable(),1.8);

    // Add the splitter to the main layout
    mainLayout->addLayout(_settingsAction.getTableSplitter());

    mainLayout->addLayout(_settingsAction.getSelectedCellClusterInfoStatusBar());
    _settingsAction.getStatusColorAction().setString("M");

    getWidget().setLayout(mainLayout);





}

void CrossSpeciesComparisonGeneDetectPlugin::adjustTableWidths(const QString& value) {
    // Assuming _settingsAction.getHorizontalLayout() returns your QHBoxLayout
    QHBoxLayout* layout = _settingsAction.getTableSplitter();
    if (!layout) return;

    QWidget* parentWidget = layout->parentWidget();
    if (!parentWidget) return;

    int totalWidth = parentWidget->width();

    double tableViewRatio = 1.2;
    double selectionDetailsTableRatio = 2;

    if (value == "small") {
        tableViewRatio = 1.8;
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
}



QColor getColorFromValue(int value, int min, int max) {
    if (value < min) value = min;
    if (value > max) value = max;

    int range = max - min;
    if (range == 0) return QColor(Qt::gray);

    int blue = 255 * (value - min) / range;

    return QColor(255 - blue, 255 - blue, 255);
}



void CrossSpeciesComparisonGeneDetectPlugin::modifyTableData()
{
    auto variant = _settingsAction.getTableModelAction().getVariant();
    QStandardItemModel* model = qobject_cast<QStandardItemModel*>(variant.value<QAbstractItemModel*>());

    if (!_settingsAction.getTableView()) {
        qDebug() << "_settingsAction.getTableView() is null";
        return;
    }

    if (!model) {
        qDebug() << "Model is null";
        QAbstractItemModel* currentModel = _settingsAction.getTableView()->model();
        if (currentModel) {
            currentModel->removeRows(0, currentModel->rowCount());
            _settingsAction.getTableView()->update();
        }
        else {
            qDebug() << "TableView model is null";
        }
        return;
    }

    _settingsAction.getTableView()->setModel(model);
    QFontMetrics metrics(_settingsAction.getTableView()->horizontalHeader()->font());
    int headerHeight = metrics.height() * 3; // Assuming 3 lines of text
    _settingsAction.getTableView()->horizontalHeader()->setFixedHeight(headerHeight);

    //QVector<int> columns = { 0,2, 3,4 };
    auto shownColumns= _settingsAction.getHiddenShowncolumns().getSelectedOptions();

    for (int i = 0; i < _settingsAction.getTableView()->model()->columnCount(); i++) {
        if (!shownColumns.contains(model->horizontalHeaderItem(i)->text())) {
            _settingsAction.getTableView()->hideColumn(i);
        }
    } 
    model->sort(1,Qt::DescendingOrder);
    //set column width
    _settingsAction.getTableView()->resizeColumnsToContents();



    //connect(_settingsAction.getTableView(), &QTableView::clicked, [this](const QModelIndex& index) {
    //    // Check if the clicked row is already selected
    //    if (_settingsAction.getTableView()->selectionModel()->isSelected(index)) {
    //        // Clear the current index if there's no selection
    //        _settingsAction.getTableView()->clearSelection();

    //        // Temporarily disable the selection mode to remove highlight
    //        QAbstractItemView::SelectionMode oldMode = _settingsAction.getTableView()->selectionMode();
    //        _settingsAction.getTableView()->setSelectionMode(QAbstractItemView::NoSelection);

    //        // Clear the current index
    //        _settingsAction.getTableView()->selectionModel()->setCurrentIndex(QModelIndex(), QItemSelectionModel::NoUpdate);

    //        // Restore the original selection mode
    //        _settingsAction.getTableView()->setSelectionMode(oldMode);
    //        // Update the view to ensure changes are reflected
    //        _settingsAction.getTableView()->update();
    //        _settingsAction.getSelctedSpeciesVals().setString("");
    //    }
    //    });



    connect(_settingsAction.getTableView()->selectionModel(), &QItemSelectionModel::currentChanged, [this](const QModelIndex& current, const QModelIndex& previous) {
        if (!current.isValid()) return;

        QString gene = current.siblingAtColumn(0).data().toString();
        _settingsAction.getSelectedGeneAction().setString(gene);
        _settingsAction.getSelectedRowIndexAction().setString(QString::number(current.row()));
        _settingsAction.getRemoveRowSelection().setEnabled(true);

        std::map<QString, Statistics> speciesExpressionMap;
        QStringList finalsettingSpeciesNamesArray;
        QString finalSpeciesNameString;
        QJsonObject valueStringReference;
        bool treeDataFound = false;
        bool isEditTreePresent = _settingsAction.getFilteringEditTreeDatasetAction().getCurrentDataset().isValid();
        const auto* model = current.model();
        const int columnCount = model->columnCount();
        const auto initColumnNamesSize = _settingsAction.getInitColumnNames().size(); 
        for (int i = 0; i < columnCount; ++i) {
            const QString columnName = model->headerData(i, Qt::Horizontal).toString();
            const auto data = current.siblingAtColumn(i).data();

            if (columnName == "Newick tree" && isEditTreePresent) {
                treeDataFound = true;
                valueStringReference = QJsonDocument::fromJson(data.toString().toUtf8()).object();
            }
            else if (columnName == "Gene Appearance Species Names") {
                finalsettingSpeciesNamesArray = data.toString().split(";");
                finalSpeciesNameString = finalsettingSpeciesNamesArray.join(" @%$,$%@ ");
            }
            else if (columnName == "Statistics") {
                // Ensure finalsettingSpeciesNamesArray is populated before processing statistics
                if (!finalsettingSpeciesNamesArray.isEmpty()) {
                    std::map<QString, Statistics> statisticsValues = convertToStatisticsMap(data.toString());
                    speciesExpressionMap = statisticsValues;
                    selectedCellCountStatusBarRemove();
                    selectedCellStatisticsStatusBarAdd(statisticsValues, finalsettingSpeciesNamesArray);
                }
                else {
                    // Handle the case where finalsettingSpeciesNamesArray is empty
                    // For example, log a warning or skip processing this column
                    qDebug() << "Warning: Gene Appearance Species Names must be processed before Statistics.";
                }
            }
        }



        std::vector<std::seed_seq::result_type> selectedSpeciesIndices;
        auto speciesDataset = _settingsAction.getSpeciesNamesDataset().getCurrentDataset();
        auto umapDataset = _settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset();
        auto mainPointsDataset = _settingsAction.getMainPointsDataset().getCurrentDataset();
        std::vector<std::seed_seq::result_type> filtSelectInndx;


        if (!speciesDataset.isValid() || !umapDataset.isValid() || !mainPointsDataset.isValid() || !_settingsAction.getFilteredUMAPDatasetPoints().isValid() || !_settingsAction.getFilteredUMAPDatasetColors().isValid())
        {
            qDebug() << "One of the datasets is not valid";
            return;
        }


        {
            auto speciesClusterDataset = mv::data().getDataset<Clusters>(speciesDataset.getDatasetId());
            auto umapPointsDataset = mv::data().getDataset<Points>(umapDataset.getDatasetId());
            auto fullMainPointsDataset = mv::data().getDataset<Points>(mainPointsDataset.getDatasetId());

            std::unordered_set<QString> speciesNamesSet(finalsettingSpeciesNamesArray.begin(), finalsettingSpeciesNamesArray.end());
            for (const auto& species : speciesClusterDataset->getClusters()) {
                if (speciesNamesSet.find(species.getName()) != speciesNamesSet.end()) {
                    const auto& indices = species.getIndices();
                    selectedSpeciesIndices.insert(selectedSpeciesIndices.end(), indices.begin(), indices.end());
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
                int tempnumPoints = selectedSpeciesIndices.size();
                int tempNumDimensions = geneIndicesSpecies.size();
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
                int tempnumPointsColors = selectedSpeciesIndices.size();
                
                std::vector<QString> columnGeneColors = { gene };
                int tempNumDimensionsColors = columnGeneColors.size();
                applyLogTransformation(resultContainerSpeciesColors);
                _settingsAction.populatePointData(speciesColorDataId, resultContainerSpeciesColors, tempnumPointsColors, tempNumDimensionsColors, columnGeneColors);

                auto scatterplotViewFactory = mv::plugins().getPluginFactory("Scatterplot View");
                mv::gui::DatasetPickerAction* colorDatasetPickerAction;
                mv::gui::DatasetPickerAction* pointDatasetPickerAction;


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
                                //auto focusSelectionAction = dynamic_cast<ToggleAction*>(plugin->findChildByPath("Focus selection"));
                                if (focusSelectionAction)
                                {
                                    focusSelectionAction->setChecked(true);

                                }
                                //"Settings/Plot/Point/Point opacity"
                                //"Settings/Plot/Point/Point opacity/Point opacity"
                                auto opacityAction = dynamic_cast<DecimalAction*>(plugin->findChildByPath("Settings/Plot/Point/Point opacity/Point opacity"));
                                if (opacityAction)
                                {
                                    opacityAction->setValue(20.0);
                                }

                            }
                        }
                    }
                }

                
               // qDebug() << "Selected species indices size: " << selectedSpeciesIndices.size();
               // qDebug()<<"datasetSize: "<<_settingsAction.getFilteredUMAPDatasetPoints()->getNumPoints();
               // qDebug()<< "filtSelectInndx points value range"<< *std::min_element(filtSelectInndx.begin(), filtSelectInndx.end()) << " " << *std::max_element(filtSelectInndx.begin(), filtSelectInndx.end());
                _settingsAction.getFilteredUMAPDatasetPoints()->setSelectionIndices(filtSelectInndx);
                mv::events().notifyDatasetDataSelectionChanged(_settingsAction.getFilteredUMAPDatasetPoints());


            }


        }









        if (treeDataFound && isEditTreePresent) {
            auto treeDataset = mv::data().getDataset<CrossSpeciesComparisonTree>(_settingsAction.getFilteringEditTreeDatasetAction().getCurrentDataset().getDatasetId());

            if (!valueStringReference.isEmpty()) {
                treeDataset->setTreeData(valueStringReference);
                events().notifyDatasetDataChanged(treeDataset);
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
            tsneDataset->setSelectionIndices(selectedPoints);
            mv::events().notifyDatasetDataSelectionChanged(tsneDataset);
        }

        if (_settingsAction.getScatterplotReembedColorOption().getCurrentText() == "Expression") {
            
            
            
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
                        indexOfGene.push_back(it - dimsValsTemp.begin());
                    }

                    std::vector<int> tempselectIndices(selectedPointsMain->getNumPoints());
                    std::iota(tempselectIndices.begin(), tempselectIndices.end(), 0);

                    if (indexOfGene.size() > 0 && tempselectIndices.size() > 0)
                    {
                        selectedPointsMain->populateDataForDimensions(resultContainerColorPoints, indexOfGene, tempselectIndices);



                        int rowSizeEmbd = rowSize;
                        int columnSizeEmbd = 1;
                        std::vector<QString> columnGeneEmbd = { gene };
                        applyLogTransformation(resultContainerColorPoints);
                        _settingsAction.populatePointData(datasetIdEmb, resultContainerColorPoints, rowSizeEmbd, columnSizeEmbd, columnGeneEmbd);

                    }
                }



            }
        }
        //_settingsAction.getSelctedSpeciesVals().setString("");
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
void CrossSpeciesComparisonGeneDetectPlugin::selectedCellCountStatusBarAdd()
{
    if (!_settingsAction.getSelectedSpeciesCellCountMap().empty())
    {
        // Create a new model for the table view
        QStandardItemModel* model = new QStandardItemModel();

        // Set headers
        model->setHorizontalHeaderLabels({ "Species", "Count\nSelected","Count\nNon\nSelected" });

        for (const auto& [species, details] : _settingsAction.getSelectedSpeciesCellCountMap()) {
            QColor backgroundColor = QColor(details.color); // Ensure color is converted to QColor

            // Calculate the brightness of the background color
            qreal brightness = backgroundColor.lightnessF();

            // Choose text color based on the brightness of the background color
            QColor textColor = (brightness > 0.5) ? Qt::black : Qt::white;

            QList<QStandardItem*> rowItems;
            QStandardItem* speciesItem = new QStandardItem(species);
            speciesItem->setBackground(backgroundColor);
            speciesItem->setForeground(textColor); // Set text color
            rowItems << speciesItem;

            QStandardItem* item;
            item = new QStandardItem();
            item->setData(QVariant(details.selectedCellsCount), Qt::EditRole);
            rowItems << item;

            item = new QStandardItem();
            item->setData(QVariant(details.nonSelectedCellsCount), Qt::EditRole);
            rowItems << item;

            model->appendRow(rowItems);
        }
        model->sort(1, Qt::DescendingOrder);
        _settingsAction.getSelectionDetailsTable()->setSelectionMode(QAbstractItemView::NoSelection);
        _settingsAction.getSelectionDetailsTable()->setModel(model);
        _settingsAction.getSelectionDetailsTable()->verticalHeader()->hide();
        _settingsAction.getSelectionDetailsTable()->resizeColumnsToContents();
        emit model->layoutChanged();

    }
    adjustTableWidths("small");
}

void CrossSpeciesComparisonGeneDetectPlugin::selectedCellStatisticsStatusBarAdd(std::map<QString, Statistics> statisticsValues, QStringList selectedSpecies)
{
    if (!_settingsAction.getSelectedSpeciesCellCountMap().empty())
    {
        // Create a new model for the table view
        QStandardItemModel* model = new QStandardItemModel();

        // Set headers
        model->setHorizontalHeaderLabels({ "Species", "Count\nSelected","Count\nNon\nSelected", "Mean\nSelected", "Mean\nNon\nSelected","Mean\nDifference","Median\nSelected", "Median\nNon\nSelected","Mode\nSelected", "Mode\nNon\nSelected","Range\nSelected", "Range\nNon\nSelected"});

        // Populate the model with sorted data and statistics
        for (const auto& [species, details] : _settingsAction.getSelectedSpeciesCellCountMap()) {
            QColor backgroundColor = QColor(details.color);

            // Calculate the brightness of the background color
            qreal brightness = backgroundColor.lightnessF();

            // Choose text color based on the brightness of the background color
            QColor textColor = (brightness > 0.5) ? Qt::black : Qt::white;

            QList<QStandardItem*> rowItems;
            QStandardItem* item = new QStandardItem(species);
            item->setBackground(backgroundColor);
            item->setForeground(textColor); // Set text color
            rowItems << item;


            // Find statistics for the species
            auto it = statisticsValues.find(species);
            if (it != statisticsValues.end()) {

                QStandardItem* item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.countSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.countNonSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.meanSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.meanNonSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                float difference = (it->second.meanSelected - it->second.meanNonSelected);
                item->setData(QVariant(QString::number(difference, 'f', 2)), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.medianSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.medianNonSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.modeSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.modeNonSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.rangeSelected), Qt::EditRole);
                rowItems << item;

                item = new QStandardItem();
                item->setData(QVariant(it->second.rangeNonSelected), Qt::EditRole);
                rowItems << item;



            }
            else {
                // Fill with placeholders if no statistics found
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
                rowItems << new QStandardItem("N/A");
            }


            // Check if the species is in the selectedSpecies list and color the row if it is
            if (selectedSpecies.contains(species)) {
                for (int i = 1; i < rowItems.size(); ++i) { // Start from 1 to skip the first column
                    rowItems[i]->setBackground(QBrush(QColor("#00A2ED")));
                }
            }



            model->appendRow(rowItems);
        }

        model->sort(5, Qt::DescendingOrder);
        _settingsAction.getSelectionDetailsTable()->setSelectionMode(QAbstractItemView::NoSelection);

        _settingsAction.getSelectionDetailsTable()->setModel(model);
        _settingsAction.getSelectionDetailsTable()->verticalHeader()->hide();
        _settingsAction.getSelectionDetailsTable()->resizeColumnsToContents();
        emit model->layoutChanged();
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


void CrossSpeciesComparisonGeneDetectPlugin::updateSpeciesData(QJsonObject& node, const std::map<QString, Statistics>& speciesExpressionMap) {
    // Check if the "name" key exists in the current node
    if (node.contains("name")) {
        QString nodeName = node["name"].toString();
        auto it = speciesExpressionMap.find(nodeName);
        // If the "name" is found in the speciesExpressionMap, update "mean" if it exists or add "mean" if it doesn't exist
        if (it != speciesExpressionMap.end()) {
            node["mean"] = std::round(it->second.meanSelected * 100.0) / 100.0; // Round to 2 decimal places

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


void CrossSpeciesComparisonGeneDetectPlugin::fromVariantMap(const QVariantMap& variantMap)
{
    ViewPlugin::fromVariantMap(variantMap);

    mv::util::variantMapMustContain(variantMap, "CSCGDV:CrossSpeciesComparison Gene Detect Plugin Settings");
    _settingsAction.fromVariantMap(variantMap["CSCGDV:CrossSpeciesComparison Gene Detect Plugin Settings"].toMap());
   // modifyTableData();
    _settingsAction.getStartComputationTriggerAction().trigger();

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

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
Q_PLUGIN_METADATA(IID "studio.manivault.CrossSpeciesComparisonGeneDetectPlugin")

using namespace mv;

CrossSpeciesComparisonGeneDetectPlugin::CrossSpeciesComparisonGeneDetectPlugin(const PluginFactory* factory) :
    ViewPlugin(factory),
    _tableView(),
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
                    if (treeDataset.isValid() && _tableView && selectedRow >= 0)
                    {
                        QString treeData = _tableView->model()->index(selectedRow, 2).data().toString();
                        //qDebug()<< "Tree data: " << treeData;
                        if (!treeData.isEmpty())
                        {

                            QJsonObject valueStringReference = QJsonDocument::fromJson(treeData.toUtf8()).object();
                            if (!valueStringReference.isEmpty())
                            {
                                treeDataset->setTreeData(valueStringReference);
                                events().notifyDatasetDataChanged(treeDataset);
                                //QString firstColumnValue = _tableView->model()->index(selectedRow, 0).data().toString();
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
                    firstColumnValues << _tableView->model()->index(row, 0).data().toString();
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
            if (_tableView && _tableView->selectionModel()) {
                // Clear the current index if there's no selection
                _tableView->clearSelection();

                // Temporarily disable the selection mode to remove highlight
                QAbstractItemView::SelectionMode oldMode = _tableView->selectionMode();
                _tableView->setSelectionMode(QAbstractItemView::NoSelection);

                // Clear the current index
                _tableView->selectionModel()->setCurrentIndex(QModelIndex(), QItemSelectionModel::NoUpdate);

                // Restore the original selection mode
                _tableView->setSelectionMode(oldMode);
                // Update the view to ensure changes are reflected
                _tableView->update();
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

        QStandardItemModel* model = qobject_cast<QStandardItemModel*>(_tableView->model());

        if (model) {
            for (int i = 0; i < model->columnCount(); i++) {
                if (!shownColumns.contains(model->horizontalHeaderItem(i)->text())) {
                    _tableView->hideColumn(i);
                }
                else
                {
                    _tableView->showColumn(i);

                }
            }
            emit model->layoutChanged();
        }
        };
    connect(&_settingsAction.getHiddenShowncolumns(), &OptionsAction::selectedOptionsChanged, this, updateHideShowColumns);



    _tableView = new QTableView();
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

    //_make the headers two three lines so that they are fully visible
    _tableView->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
    _tableView->horizontalHeader()->setStretchLastSection(true);
    _tableView->horizontalHeader()->setMinimumSectionSize(50);
    _tableView->horizontalHeader()->setMaximumSectionSize(600);
    _tableView->horizontalHeader()->setHighlightSections(false);
    _tableView->horizontalHeader()->setSortIndicatorShown(true);
    //change height of headers

    

    //make long strings in the cells visible and not ...shortened
    //_tableView->setTextElideMode(Qt::ElideNone);
    //_tableView->setWordWrap(true);
    //_tableView->setAlternatingRowColors(true);
    //_tableView->setSortingEnabled(true);

    //on hovering a cell, show the full text available in a tooltip
    connect(_tableView, &QTableView::entered, [this](const QModelIndex& index) {
        if (index.isValid()) {
            QString text = index.data().toString();
            if (!text.isEmpty()) {
                _tableView->setToolTip(text);
            }
        }
        });


    /*
    connect(_tableView, &QTableView::clicked, [this](const QModelIndex& index) {
        QModelIndex firstColumnIndex = index.sibling(index.row(), 0);
        auto gene = firstColumnIndex.data().toString();
        _settingsAction.getSelectedGeneAction().setString(gene);

        //if (QApplication::keyboardModifiers() & Qt::ShiftModifier) 
        
        //{
            // If Shift is pressed, add the row to the selection
          //  _tableView->selectionModel()->select(index, QItemSelectionModel::Select | QItemSelectionModel::Rows);
        //}
        //else {
            // If Shift is not pressed, select only this row
            //_tableView->selectionModel()->clearSelection();
           // _tableView->selectionModel()->select(index, QItemSelectionModel::Select | QItemSelectionModel::Rows);
       // }

        // Get the selected rows and convert them to a string list
        //QModelIndexList selectedRows = _tableView->selectionModel()->selectedRows();
        QStringList selectedRowsStrList;
        for (const QModelIndex& selectedIndex : selectedRows) {
            selectedRowsStrList << QString::number(selectedIndex.row());
        }

        // Join the string list into a single string with comma separation
        QString selectedRowsStr = selectedRowsStrList.join(",");
        _settingsAction.getSelectedRowIndexAction().setString(selectedRowsStr);
        });
   */


    _tableView->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _tableView->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _tableView->sortByColumn(3, Qt::DescendingOrder);

    _tableView->verticalHeader()->hide();
    _tableView->setMouseTracking(true);
    _tableView->setToolTipDuration(10000);
    QFont font = _tableView->horizontalHeader()->font();
    font.setBold(true);
    _tableView->horizontalHeader()->setFont(font);
    _tableView->setStyleSheet("QTableView::item:selected { background-color: #00A2ED; }");
    _tableView->horizontalHeader()->setHighlightSections(false);
    _tableView->verticalHeader()->setHighlightSections(false);


    auto mainLayout = new QVBoxLayout();
    mainLayout->setContentsMargins(0, 0, 0, 0);
    mainLayout->setSpacing(0);

    auto mainOptionsLayout = new QHBoxLayout();
    mainOptionsLayout->setSpacing(0);
    mainOptionsLayout->setContentsMargins(0, 0, 0, 0);
    auto extraOptionsGroup= new VerticalGroupAction(this,"Settings");

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
    
    auto tsneOptionsGroup= new VerticalGroupAction(this,"Options");
    tsneOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("tools"));
    tsneOptionsGroup->addAction(&_settingsAction.getUsePreComputedTSNE());
    tsneOptionsGroup->addAction(&_settingsAction.getTsnePerplexity());
    tsneOptionsGroup->addAction(&_settingsAction.getHiddenShowncolumns());

    auto mainOptionsGroupLayout = new QVBoxLayout();
    auto mainOptionsGroup1 = new HorizontalGroupAction(this, "MainGroup1");
    auto mainOptionsGroup2 = new HorizontalGroupAction(this, "MainGroup2");
    mainOptionsGroup1->setIcon(Application::getIconFont("FontAwesome").getIcon("database"));
    mainOptionsGroup2->setIcon(Application::getIconFont("FontAwesome").getIcon("play"));
    mainOptionsGroup2->addAction(&_settingsAction.getStartComputationTriggerAction());
    mainOptionsGroup2->addAction(&_settingsAction.getRemoveRowSelection());
    mainOptionsGroup2->addAction(&_settingsAction.getScatterplotReembedColorOption());
    mainOptionsGroup1->addAction(&_settingsAction.getTopNGenesFilter());
    mainOptionsGroup1->addAction(&_settingsAction.getTypeofTopNGenes());
    auto group1Widget= mainOptionsGroup1->createWidget(&getWidget());
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

    //fullSettingsLayout->addWidget(_settingsAction.getSelectedCellClusterInfoStatusBar());

    mainLayout->addLayout(fullSettingsLayout);
    mainLayout->addWidget(_tableView);
    _settingsAction.getSelectedCellClusterInfoStatusBar()->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    mainLayout->addWidget(_settingsAction.getSelectedCellClusterInfoStatusBar());
    _settingsAction.getStatusColorAction().setString("M");
    // Set the layout for the widget
    getWidget().setLayout(mainLayout);





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
    QStandardItemModel* model = variant.value<QStandardItemModel*>();

    if (_tableView == nullptr) {
        qDebug() << "_tableView is null";
        return;
    }

    if (model == nullptr) {
        qDebug() << "Model is null";
        if (_tableView->model() != nullptr) {
            _tableView->model()->removeRows(0, _tableView->model()->rowCount());
            _tableView->update();
        }
        else {
            qDebug() << "TableView model is null";
        }
        return;
    }

    _tableView->setModel(model);


    //QVector<int> columns = { 0,2, 3,4 };
    auto shownColumns= _settingsAction.getHiddenShowncolumns().getSelectedOptions();

    for (int i = 0; i < _tableView->model()->columnCount(); i++) {
        if (!shownColumns.contains(model->horizontalHeaderItem(i)->text())) {
            _tableView->hideColumn(i);
        }
    } 
    model->sort(3,Qt::DescendingOrder);


    //connect(_tableView, &QTableView::clicked, [this](const QModelIndex& index) {
    //    // Check if the clicked row is already selected
    //    if (_tableView->selectionModel()->isSelected(index)) {
    //        // Clear the current index if there's no selection
    //        _tableView->clearSelection();

    //        // Temporarily disable the selection mode to remove highlight
    //        QAbstractItemView::SelectionMode oldMode = _tableView->selectionMode();
    //        _tableView->setSelectionMode(QAbstractItemView::NoSelection);

    //        // Clear the current index
    //        _tableView->selectionModel()->setCurrentIndex(QModelIndex(), QItemSelectionModel::NoUpdate);

    //        // Restore the original selection mode
    //        _tableView->setSelectionMode(oldMode);
    //        // Update the view to ensure changes are reflected
    //        _tableView->update();
    //        _settingsAction.getSelctedSpeciesVals().setString("");
    //    }
    //    });



    connect(_tableView->selectionModel(), &QItemSelectionModel::currentChanged, [this](const QModelIndex& current, const QModelIndex& previous) {
        if (!current.isValid()) return;

        QString gene = current.siblingAtColumn(0).data().toString();
        _settingsAction.getSelectedGeneAction().setString(gene);
        _settingsAction.getSelectedRowIndexAction().setString(QString::number(current.row()));
        _settingsAction.getRemoveRowSelection().setEnabled(true);

        std::map<QString, float> speciesExpressionMap;
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
            //for all columns other than _settingsAction.getInitColumnNames().size(
            if (i > initColumnNamesSize) {
                speciesExpressionMap[columnName] = data.toFloat();
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
                        _settingsAction.populatePointData(datasetIdEmb, resultContainerColorPoints, rowSizeEmbd, columnSizeEmbd, columnGeneEmbd);

                    }
                }



            }
        }
        _settingsAction.getSelctedSpeciesVals().setString(finalSpeciesNameString);
        });



    emit model->layoutChanged();

}
void CrossSpeciesComparisonGeneDetectPlugin::updateSpeciesData(QJsonObject& node, const std::map<QString, float>& speciesExpressionMap) {
    // Check if the "name" key exists in the current node
    if (node.contains("name")) {
        QString nodeName = node["name"].toString();
        auto it = speciesExpressionMap.find(nodeName);
        // If the "name" is found in the speciesExpressionMap, update "mean" if it exists or add "mean" if it doesn't exist
        if (it != speciesExpressionMap.end()) {
            node["mean"] = it->second; // Use it->second to access the value in the map
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

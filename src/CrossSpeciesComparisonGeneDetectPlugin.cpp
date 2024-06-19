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

            if (_settingsAction.getFilteringTreeDatasetAction().getCurrentDataset().isValid())
            {
                auto treeDataset = mv::data().getDataset<CrossSpeciesComparisonTree>(_settingsAction.getFilteringTreeDatasetAction().getCurrentDataset().getDatasetId());
              
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
            else
            {
                qDebug() << "Tree dataset is not valid";
            }
        };

    connect(&_settingsAction.getSelectedRowIndexAction(), &StringAction::stringChanged, this, updateSelectedRowIndex);

    const auto updateSelectedGene = [this]() -> void
        {


        };

    connect(&_settingsAction.getSelectedGeneAction(), &StringAction::stringChanged, this, updateSelectedGene);

    const auto removeRowSelectionTable = [this]() -> void
        {
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


                if (_settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset().isValid() && _settingsAction.getScatterplotEmbeddingColorOption().getCurrentDataset().isValid())
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
                                        colorDatasetPickerAction->setCurrentDataset(_settingsAction.getScatterplotEmbeddingColorOption().getCurrentDataset());

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
    extraOptionsGroup->addAction(&_settingsAction.getFilteringTreeDatasetAction());
    extraOptionsGroup->addAction(&_settingsAction.getOptionSelectionAction());
    extraOptionsGroup->addAction(&_settingsAction.getReferenceTreeDatasetAction());
    extraOptionsGroup->addAction(&_settingsAction.getMainPointsDataset());
    extraOptionsGroup->addAction(&_settingsAction.getEmbeddingDataset());
    extraOptionsGroup->addAction(&_settingsAction.getSpeciesNamesDataset());
    extraOptionsGroup->addAction(&_settingsAction.getClusterNamesDataset());
    extraOptionsGroup->addAction(&_settingsAction.getFilteredGeneNames());
    extraOptionsGroup->addAction(&_settingsAction.getGeneNamesConnection());
    extraOptionsGroup->addAction(&_settingsAction.getTsnePerplexity());
    extraOptionsGroup->addAction(&_settingsAction.getCreateRowMultiSelectTree());
    extraOptionsGroup->addAction(&_settingsAction.getPerformGeneTableTsneAction());
    extraOptionsGroup->addAction(&_settingsAction.getHiddenShowncolumns());
    extraOptionsGroup->addAction(&_settingsAction.getSelctedSpeciesVals());
    extraOptionsGroup->addAction(&_settingsAction.getScatterplotEmbeddingColorOption());
    extraOptionsGroup->addAction(&_settingsAction.getScatterplotEmbeddingPointsUMAPOption());
    


    auto mainOptionsGroupLayout = new QVBoxLayout();
    auto mainOptionsGroup1 = new HorizontalGroupAction(this, "MainGroup1");
    auto mainOptionsGroup2 = new HorizontalGroupAction(this, "MainGroup2");

    mainOptionsGroup1->setIcon(Application::getIconFont("FontAwesome").getIcon("database"));
    mainOptionsGroup2->setIcon(Application::getIconFont("FontAwesome").getIcon("play"));


    mainOptionsGroup2->addAction(&_settingsAction.getStatusAction(), -1, [this](WidgetAction* action, QWidget* widget) -> void
        {
            auto labelWidget = widget->findChild<QLabel*>("Label");

            if (labelWidget)
            {
                // Set initial state text and color
                labelWidget->setText("");
                labelWidget->setStyleSheet("background-color: none; color: white;"); // Set initial text color to white
                qDebug() << "Initial status color: " << _settingsAction.getStatusColorAction().getString();
                connect(&_settingsAction.getStatusColorAction(), &StringAction::stringChanged, this, [this, labelWidget](const QString& string) -> void
                    {
                        qDebug() << "Status color changed to: " << string;
                        QString labelText = "";
                        QString backgroundColor = "none";
                        if (string == "C")
                        {
                            labelText = "Up-to-date";
                            backgroundColor = "#28a745";
                        }
                        else if (string == "M")
                        {
                            labelText = "Outdated";
                            backgroundColor = "#ffc107";
                        }
                        else
                        {
                            labelText = "Unknown";
                            backgroundColor = "#6c757d";
                        }
                        labelWidget->setText(labelText);
                        labelWidget->setStyleSheet(QString("background-color: %1; color: white;").arg(backgroundColor));
                    });
            }

        });

    mainOptionsGroup2->addAction(&_settingsAction.getStartComputationTriggerAction());

    mainOptionsGroup2->addAction(&_settingsAction.getRemoveRowSelection());

    mainOptionsGroup1->addAction(&_settingsAction.getTopNGenesFilter());
    mainOptionsGroup1->addAction(&_settingsAction.getScatterplotReembedColorOption());

    auto group1Widget= mainOptionsGroup1->createWidget(&getWidget());
    group1Widget->setMaximumWidth(550);
    mainOptionsGroupLayout->addWidget(group1Widget);

    auto group2Widget = mainOptionsGroup2->createWidget(&getWidget());
    group2Widget->setMaximumWidth(400);
    mainOptionsGroupLayout->addWidget(group2Widget);  

    mainOptionsLayout->addLayout(mainOptionsGroupLayout);
    mainOptionsLayout->addWidget(extraOptionsGroup->createCollapsedWidget(&getWidget()), 1);
    
    mainLayout->addLayout(mainOptionsLayout);



    //
    if (0)
    {
        // Create a new QSplitter
        QSplitter* splitter = new QSplitter();

        // Add _tableView to the splitter
        splitter->addWidget(_tableView);

        // Create another view
        QWidget* anotherView = new QWidget();
        splitter->addWidget(anotherView);

        // Set stretch factors for the widgets
        splitter->setStretchFactor(0, 1); // _tableView
        splitter->setStretchFactor(1, 1); // anotherView

        // Get the total available width
        int totalWidth = splitter->width();

        // Calculate the width for each widget
        int widgetWidth = totalWidth / 2; // divide by the number of widgets

        // Set the sizes of the child widgets
        QList<int> sizes;
        sizes << widgetWidth << widgetWidth; // adjust these values as needed
        splitter->setSizes(sizes);

        // Set the splitter as the main widget in your layout
        mainLayout->addWidget(splitter);
    }
    else
    {
        mainLayout->addWidget(_tableView);
    }


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
        




        std::map<QString, float> speciesExpressionMap;
        QStringList finalsettingSpeciesNamesArray;
        QString finalSpeciesNameString;
        QJsonObject valueStringReference;
        bool treeDataFound = false;

        const auto* model = current.model();
        const int columnCount = model->columnCount();
        for (int i = 0; i < columnCount; ++i) {
            const QString columnName = model->headerData(i, Qt::Horizontal).toString();
            const auto data = current.siblingAtColumn(i).data();
            if (i > 5) {
                speciesExpressionMap[columnName] = data.toFloat();
            }
            else if (columnName == "Newick tree") {
                treeDataFound = true;
                valueStringReference = QJsonDocument::fromJson(data.toString().toUtf8()).object();
            }
            else if (columnName == "Gene Apearance Species Names") {
                finalsettingSpeciesNamesArray = data.toString().split(";");
                finalSpeciesNameString = finalsettingSpeciesNamesArray.join(" @%$,$%@ ");
            }
        }



        std::vector<std::seed_seq::result_type> selectedSpeciesIndices;
        auto speciesDataset = _settingsAction.getSpeciesNamesDataset().getCurrentDataset();
        auto umapDataset = _settingsAction.getScatterplotEmbeddingPointsUMAPOption().getCurrentDataset();
        auto mainPointsDataset = _settingsAction.getMainPointsDataset().getCurrentDataset();
        
        if (speciesDataset.isValid() && umapDataset.isValid() && mainPointsDataset.isValid() && _settingsAction.getFilteredUMAPDatasetPoints().isValid() && _settingsAction.getFilteredUMAPDatasetColors().isValid())
        {
            auto speciesClusterDataset = mv::data().getDataset<Clusters>(speciesDataset.getDatasetId());
            auto umapPointsDataset = mv::data().getDataset<Points>(umapDataset.getDatasetId());
            auto fullMainPointsDataset = mv::data().getDataset<Points>(mainPointsDataset.getDatasetId());

            auto speciesClusters= speciesClusterDataset->getClusters();
            for (const auto& species : speciesClusters) {
                if (finalsettingSpeciesNamesArray.contains(species.getName())) {
                    const auto& indices = species.getIndices();
                    selectedSpeciesIndices.insert(selectedSpeciesIndices.end(), indices.begin(), indices.end());
                }
            }
            //
           //getFilteredUMAPDatasetColors()
            auto dimensionNamesUmap = umapPointsDataset->getDimensionNames();
            std::vector<int> geneIndicesSpecies;
            for (int i = 0; i < umapPointsDataset->getNumDimensions(); i++)
            {
                geneIndicesSpecies.push_back(i);
            }

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

                for (int i = 0; i < fullMainPointsDataset->getDimensionNames().size(); i++)
                {
                    if (fullMainPointsDataset->getDimensionNames()[i] == gene)
                    {
                        selectedGeneIndex.push_back(i);
                        break;
                    }
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
                            }
                        }
                    }
                }
               // _settingsAction.getFilteredUMAPDatasetPoints()->setSelectionIndices(_settingsAction.getSelectedIndicesFromStorage());
                //mv::events().notifyDatasetDataSelectionChanged(_settingsAction.getFilteredUMAPDatasetPoints());


            }


        }









        if (treeDataFound) {
            auto treeDataset = mv::data().getDataset<CrossSpeciesComparisonTree>(_settingsAction.getFilteringTreeDatasetAction().getCurrentDataset().getDatasetId());
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
            if (speciesColorClusterDataset.isValid() && expressionColorPointDataset.isValid()) {
                const int rowSize = expressionColorPointDataset->getNumPoints();
                std::vector<float> resultContainerColorPoints(rowSize, -1.0);
                const QString datasetIdEmb = expressionColorPointDataset->getId();

                for (const auto& species : speciesColorClusterDataset->getClusters()) {
                    float speciesValue = speciesExpressionMap[species.getName()];
                    for (auto index : species.getIndices()) {
                        resultContainerColorPoints[index] = speciesValue;
                    }
                }

                QString tempDatasetIdEmb = datasetIdEmb; // Assuming datasetIdEmb is const QString
                int rowSizeEmbd = rowSize;
                int columnSizeEmbd = 1;
                std::vector<QString> columnGeneEmbd = { gene };
                _settingsAction.populatePointData(tempDatasetIdEmb, resultContainerColorPoints, rowSizeEmbd, columnSizeEmbd, columnGeneEmbd);

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

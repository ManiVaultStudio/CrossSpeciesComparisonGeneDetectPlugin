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
                _settingsAction.getGeneNamesConnection().setString(firstColumnValue);


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

    const auto updateTableModel = [this]() -> void
        {
            modifyTableData();

        };

    connect(&_settingsAction.getTableModelAction(), &VariantAction::variantChanged, this, updateTableModel);

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
    _tableView->setSelectionMode(QAbstractItemView::ExtendedSelection);
    //only highlight multiple rows if shiuft is pressed
    _tableView->setSelectionBehavior(QAbstractItemView::SelectRows);

    connect(_tableView, &QTableView::clicked, [this](const QModelIndex& index) {
        QModelIndex firstColumnIndex = index.sibling(index.row(), 0);
        auto gene = firstColumnIndex.data().toString();
        _settingsAction.getSelectedGeneAction().setString(gene);

        if (QApplication::keyboardModifiers() & Qt::ShiftModifier) {
            // If Shift is pressed, add the row to the selection
            _tableView->selectionModel()->select(index, QItemSelectionModel::Select | QItemSelectionModel::Rows);
        }
        else {
            // If Shift is not pressed, select only this row
            _tableView->selectionModel()->clearSelection();
            _tableView->selectionModel()->select(index, QItemSelectionModel::Select | QItemSelectionModel::Rows);
        }

        // Get the selected rows and convert them to a string list
        QModelIndexList selectedRows = _tableView->selectionModel()->selectedRows();
        QStringList selectedRowsStrList;
        for (const QModelIndex& selectedIndex : selectedRows) {
            selectedRowsStrList << QString::number(selectedIndex.row());
        }

        // Join the string list into a single string with comma separation
        QString selectedRowsStr = selectedRowsStrList.join(",");
        _settingsAction.getSelectedRowIndexAction().setString(selectedRowsStr);
        });
    //add a lambda function 


    _tableView->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _tableView->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _tableView->sortByColumn(1, Qt::DescendingOrder);
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
    extraOptionsGroup->addAction(&_settingsAction.getSpeciesNamesDataset());
    extraOptionsGroup->addAction(&_settingsAction.getFilteredGeneNames());
    extraOptionsGroup->addAction(&_settingsAction.getGeneNamesConnection());
    extraOptionsGroup->addAction(&_settingsAction.getTsnePerplexity());


    


    auto mainOptionsGroup = new HorizontalGroupAction(this, "Trigger");
    mainOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("play"));
    mainOptionsGroup->addAction(&_settingsAction.getStartComputationTriggerAction());
    mainOptionsGroup->addAction(&_settingsAction.getTopNGenesFilter());
    mainOptionsGroup->addAction(&_settingsAction.getCreateRowMultiSelectTree());
    mainOptionsGroup->addAction(&_settingsAction.getPerformGeneTableTsneAction());

    mainOptionsLayout->addWidget(mainOptionsGroup->createWidget(&getWidget()),2);
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
    _tableView->sortByColumn(1, Qt::DescendingOrder);

    QVector<int> columns = { 0, 1, 3,4,5,6,7 };
    for (int i = 0; i < _tableView->model()->columnCount(); i++) {
        if (!columns.contains(i)) {
            _tableView->hideColumn(i);
        }
    }
    emit model->layoutChanged();

    std::vector<float> items;
    std::map<QString, std::vector<std::seed_seq::result_type>> clusterMap;
    items.reserve(2 * model->rowCount()); // Reserve space for items

    for (int i = 0; i < model->rowCount(); i++) {
        auto index1 = model->index(i, 1);
        auto index4 = model->index(i, 4);
        auto index2 = model->index(i, 3);

        if (index1.isValid() && index4.isValid()) {
            bool ok1, ok4;
            float item1 = index1.data().toFloat(&ok1);
            float item4 = index4.data().toFloat(&ok4);
            item4= item4 * 100;// change this to a percentage
            if (ok1 && ok4) {
                items.push_back(item1);
                items.push_back(item4);
                //qDebug() << "Item1: " << item1 << " Item4: " << item4;
            }
        }

        if (index2.isValid()) {
            QString gene = index2.data().toString();
            if (!gene.isEmpty()) {
                clusterMap[gene].push_back(i);
            }
        }
    }
    //iterate clustermap and get the min and max values
    int minVal = INT_MAX;
    int maxVal = INT_MIN;
    for (const auto& [gene, indices] : clusterMap) {
        int numberVal = gene.split("/")[0].toInt();
        if (numberVal) {
            if (numberVal < minVal) minVal = numberVal;
            if (numberVal > maxVal) maxVal = numberVal;
        }
    }

    std::vector<QString> dimensionNames = { "Variance","Similarity" };
    if (!_pointsDataset.isValid() || !_clusterDataset.isValid())
    {
        
        _pointsDataset = mv::data().createDataset("Points", "GeneSimilarityDataset");
        _pointsDataset->setGroupIndex(8);
        mv::events().notifyDatasetAdded(_pointsDataset);
        _clusterDataset = mv::data().createDataset("Cluster", "GeneSimilarityClusterDataset", _pointsDataset);
        _clusterDataset->setGroupIndex(8);
        mv::events().notifyDatasetAdded(_clusterDataset);

   }

    

    if (_pointsDataset.isValid() && _clusterDataset.isValid())
    {

        _pointsDataset->setData(items.data(), items.size() / 2, 2);
        _pointsDataset->setDimensionNames(dimensionNames);
        bool canPerformTSNE = false;
        auto perplexityValue = _settingsAction.getTsnePerplexity().getValue();
        if ((items.size() / 2)> perplexityValue && perplexityValue>0)
        {
            canPerformTSNE = true;
        }
        events().notifyDatasetDataChanged(_pointsDataset);

        _clusterDataset->getClusters() = QVector<Cluster>();
        events().notifyDatasetDataChanged(_clusterDataset);
        for (const auto& [gene, indices] : clusterMap) {
            Cluster cluster;
            cluster.setIndices(indices);
            cluster.setName(gene);
            int numberVal = gene.split("/")[0].toInt();
            QColor color;

            if (numberVal) {
                color = getColorFromValue(numberVal, minVal, maxVal);
            }
            else {
                color = QColor(Qt::gray);
            }

            cluster.setColor(color);

 
            


            _clusterDataset->getClusters().append(cluster);
        }
        events().notifyDatasetDataChanged(_clusterDataset);

        if (_settingsAction.getPerformGeneTableTsneAction().isChecked() && canPerformTSNE)
        {
            if (model->columnCount()>10)
            {
                std::vector<float> allitems;

                allitems.reserve((model->columnCount() - 8)* model->rowCount()); // Reserve space for items

                for (int i = 0; i < model->rowCount(); i++)
                {
                    for (int j = 8; j < (model->columnCount() - 8); j++)
                    {
                        auto index = model->index(i, j);
                        float item = index.data().toFloat();
                        allitems.push_back(item);
                    }
                }



                _pointsDataset->setData(allitems.data(), items.size() / 2, 2);
                _pointsDataset->setDimensionNames(dimensionNames);

                events().notifyDatasetDataChanged(_pointsDataset);
            }

        
        auto analysisPlugin = mv::plugins().requestPlugin<AnalysisPlugin>("tSNE Analysis", { _pointsDataset });
        if (!analysisPlugin) {
            qDebug() << "Could not find create TSNE Analysis";
            return;
        }

        if (_lowDimTSNEDataset.isValid())
        {
            auto runningAction = dynamic_cast<TriggerAction*>(_lowDimTSNEDataset->findChildByPath("TSNE/TsneComputationAction/Running"));

            if (runningAction)
            {
                auto perplexityAction= dynamic_cast<IntegralAction*>(_lowDimTSNEDataset->findChildByPath("TSNE/Perplexity"));
                if (perplexityAction)
                {
                    qDebug()<< "Perplexity: Found" ;
                    perplexityAction->setValue(perplexityValue);
                }
                else
                {
                    qDebug() << "Perplexity: Not Found";
                }
                if (runningAction->isChecked())
                {
                    auto stopAction = dynamic_cast<TriggerAction*>(_lowDimTSNEDataset->findChildByPath("TSNE/TsneComputationAction/Stop"));
                    if (stopAction)
                    {
                        stopAction->trigger();
                        std::this_thread::sleep_for(std::chrono::seconds(5));
                    }
                }

            }
        }

        if (_lowDimTSNEDataset.isValid())
        {
            auto datasetIDLowRem = _lowDimTSNEDataset.getDatasetId();
            mv::events().notifyDatasetAboutToBeRemoved(_lowDimTSNEDataset);
            mv::data().removeDataset(_lowDimTSNEDataset);
            mv::events().notifyDatasetRemoved(datasetIDLowRem, PointType);
        }
        _lowDimTSNEDataset = analysisPlugin->getOutputDataset();
        if (_lowDimTSNEDataset.isValid())
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

                            pointDatasetPickerAction->setCurrentDataset(_lowDimTSNEDataset);

                            colorDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Color"));
                            if (colorDatasetPickerAction)
                            {
                                colorDatasetPickerAction->setCurrentText("");
                                colorDatasetPickerAction->setCurrentDataset(_clusterDataset);
                            }
                        }
                    }
                }
            }

            auto startAction = dynamic_cast<TriggerAction*>(_lowDimTSNEDataset->findChildByPath("TSNE/TsneComputationAction/Start"));
            if (startAction) {

                startAction->trigger();

                analysisPlugin->getOutputDataset()->setSelectionIndices({});
            }

        }

    }
        else
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

                            pointDatasetPickerAction->setCurrentDataset(_pointsDataset);

                            colorDatasetPickerAction = dynamic_cast<DatasetPickerAction*>(plugin->findChildByPath("Settings/Datasets/Color"));
                            if (colorDatasetPickerAction)
                            {
                                colorDatasetPickerAction->setCurrentText("");
                                colorDatasetPickerAction->setCurrentDataset(_clusterDataset);
                            }
                        }
                    }
                }
            }
        }
    
    }
    else
    {
        qDebug() << "Low dimensional dataset is not valid";

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

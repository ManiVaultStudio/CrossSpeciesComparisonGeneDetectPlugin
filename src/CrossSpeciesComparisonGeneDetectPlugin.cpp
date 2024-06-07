#include "CrossSpeciesComparisonGeneDetectPlugin.h"

#include <event/Event.h>
#include <CrossSpeciesComparisonTreeData.h>
#include <DatasetsMimeData.h>
#include <QHeaderView> 
#include <QDebug>
#include <QMimeData>
#include <QShortcut>

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
              
               int selectedRow= _settingsAction.getSelectedRowIndexAction().getString().toInt();

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
                       }
                   }
               }


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
    
    connect(_tableView, &QTableView::clicked, [this](const QModelIndex& index) {
        QModelIndex firstColumnIndex = index.sibling(index.row(), 0);
        auto gene = firstColumnIndex.data().toString();
       _settingsAction.getSelectedGeneAction().setString(gene);
       _settingsAction.getSelectedRowIndexAction().setString(QString::number(index.row()));
        });
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
    extraOptionsGroup->addAction(&_settingsAction.getTopNGenesFilter());
    

    


    auto mainOptionsGroup = new HorizontalGroupAction(this, "Trigger");
    mainOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("play"));
    mainOptionsGroup->addAction(&_settingsAction.getStartComputationTriggerAction());

    mainOptionsLayout->addWidget(mainOptionsGroup->createWidget(&getWidget()),2);
    mainOptionsLayout->addWidget(extraOptionsGroup->createCollapsedWidget(&getWidget()), 1);
    
    mainLayout->addLayout(mainOptionsLayout);
    mainLayout->addWidget(_tableView);

    getWidget().setLayout(mainLayout);



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

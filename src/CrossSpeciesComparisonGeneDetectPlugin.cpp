#include "CrossSpeciesComparisonGeneDetectPlugin.h"

#include <event/Event.h>

#include <DatasetsMimeData.h>
#include <QHeaderView> /
#include <QDebug>
#include <QMimeData>

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
    auto layout = new QVBoxLayout();

    layout->setContentsMargins(0, 0, 0, 0);
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

    //show a thin x and y axis scrollbar
    _tableView->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _tableView->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    _tableView->sortByColumn(1, Qt::DescendingOrder);
    //hide _tableView now headers
    //_tableView->horizontalHeader()->hide();
    _tableView->verticalHeader()->hide();

    //add row selection color
    _tableView->setStyleSheet("QTableView::item:selected { background-color: #00A2ED; }");
    //do not alternate row color
    //_tableView->setAlternatingRowColors(false);

    //do not highlight header
    _tableView->horizontalHeader()->setHighlightSections(false);
    _tableView->verticalHeader()->setHighlightSections(false);

    QWidget* widget = new QWidget();

    layout->addWidget(_settingsAction.getOptionSelectionAction().createWidget(&getWidget()));
    layout->addWidget(_tableView);

    getWidget().setLayout(layout);






    //QStandardItemModel* model = new QStandardItemModel(4, 4, this);
    //QStandardItemModel* model = variant.value<QStandardItemModel*>();
    // 
    // 
    ////add header 
    //model->setHorizontalHeaderItem(0, new QStandardItem("Gene"));
    //model->setHorizontalHeaderItem(1, new QStandardItem("Variance"));
    //int numOfSpecies = 25;
    //for (int i=0+2;i< numOfSpecies+2;i++)
    //{
    //    model->setHorizontalHeaderItem(i, new QStandardItem(QString("Mean_Species") + QString::number(i)));
    //}

    ////add dummy data
    ////model->setItem(0, 0, new QStandardItem("Gene1"));
    ////model->setItem(0, 1, new QStandardItem("0.5"));
    ////model->setItem(0, 2, new QStandardItem("0.1"));
    ////model->setItem(0, 3, new QStandardItem("0.2"));

    ////model->setItem(1, 0, new QStandardItem("Gene2"));
    ////model->setItem(1, 1, new QStandardItem("0.6"));
    ////model->setItem(1, 2, new QStandardItem("0.2"));
    ////model->setItem(1, 3, new QStandardItem("0.3"));

    ////model->setItem(2, 0, new QStandardItem("Gene3"));
    ////model->setItem(2, 1, new QStandardItem("0.7"));
    ////model->setItem(2, 2, new QStandardItem("0.3"));
    ////model->setItem(2, 3, new QStandardItem("0.4"));

    ////model->setItem(3, 0, new QStandardItem("Gene4"));
    ////model->setItem(3, 1, new QStandardItem("0.8"));
    ////model->setItem(3, 2, new QStandardItem("0.4"));
    ////model->setItem(3, 3, new QStandardItem("0.5"));


   
    //layout->addWidget(_currentDatasetNameLabel);

    // Apply the layout
    

    //// Instantiate new drop widget
    //_dropWidget = new DropWidget(_currentDatasetNameLabel);

    //// Set the drop indicator widget (the widget that indicates that the view is eligible for data dropping)
    //_dropWidget->setDropIndicatorWidget(new DropWidget::DropIndicatorWidget(&getWidget(), "No data loaded", "Drag an item from the data hierarchy and drop it here to visualize data..."));

    //// Initialize the drop regions
    //_dropWidget->initialize([this](const QMimeData* mimeData) -> DropWidget::DropRegions {
    //    // A drop widget can contain zero or more drop regions
    //    DropWidget::DropRegions dropRegions;

    //    const auto datasetsMimeData = dynamic_cast<const DatasetsMimeData*>(mimeData);

    //    if (datasetsMimeData == nullptr)
    //        return dropRegions;

    //    if (datasetsMimeData->getDatasets().count() > 1)
    //        return dropRegions;

    //    // Gather information to generate appropriate drop regions
    //    const auto dataset = datasetsMimeData->getDatasets().first();
    //    const auto datasetGuiName = dataset->getGuiName();
    //    const auto datasetId = dataset->getId();
    //    const auto dataType = dataset->getDataType();
    //    const auto dataTypes = DataTypes({ PointType });

    //    // Visually indicate if the dataset is of the wrong data type and thus cannot be dropped
    //    if (!dataTypes.contains(dataType)) {
    //        dropRegions << new DropWidget::DropRegion(this, "Incompatible data", "This type of data is not supported", "exclamation-circle", false);
    //    }
    //    else {

    //        // Get points dataset from the core
    //        auto candidateDataset = mv::data().getDataset<Points>(datasetId);

    //        // Accept points datasets drag and drop
    //        if (dataType == PointType) {
    //            const auto description = QString("Load %1 into example view").arg(datasetGuiName);

    //            if (_points == candidateDataset) {

    //                // Dataset cannot be dropped because it is already loaded
    //                dropRegions << new DropWidget::DropRegion(this, "Warning", "Data already loaded", "exclamation-circle", false);
    //            }
    //            else {

    //                // Dataset can be dropped
    //                dropRegions << new DropWidget::DropRegion(this, "Points", description, "map-marker-alt", true, [this, candidateDataset]() {
    //                    _points = candidateDataset;
    //                });
    //            }
    //        }
    //    }

    //    return dropRegions;
    //});

    //// Respond when the name of the dataset in the dataset reference changes
    //connect(&_points, &Dataset<Points>::guiNameChanged, this, [this]() {

    //    auto newDatasetName = _points->getGuiName();

    //    // Update the current dataset name label
    //    _currentDatasetNameLabel->setText(QString("Current points dataset: %1").arg(newDatasetName));

    //    // Only show the drop indicator when nothing is loaded in the dataset reference
    //    _dropWidget->setShowDropIndicator(newDatasetName.isEmpty());
    //});

    //// Alternatively, classes which derive from hdsp::EventListener (all plugins do) can also respond to events
    //_eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DatasetAdded));
    //_eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DatasetDataChanged));
    //_eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DatasetRemoved));
    //_eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DatasetDataSelectionChanged));
    //_eventListener.registerDataEventByType(PointType, std::bind(&CrossSpeciesComparisonGeneDetectPlugin::onDataEvent, this, std::placeholders::_1));
}

void CrossSpeciesComparisonGeneDetectPlugin::modifyTableData()
{
    auto variant = _settingsAction.getTableModelAction().getVariant();
    // variant to QStandardItemModel
    QStandardItemModel* model = variant.value<QStandardItemModel*>();
    if (_tableView != nullptr) {
        if (model != nullptr) {
            _tableView->setModel(model);
        }
        else {
            // Handle the case where model is null
            qDebug() << "Model is null";
            if (_tableView->model() != nullptr) {
                _tableView->model()->removeRows(0, _tableView->model()->rowCount());
            }
            else {
                qDebug() << "TableView model is null";
            }
        }
    }
    else {
        qDebug() << "_tableView is null";
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

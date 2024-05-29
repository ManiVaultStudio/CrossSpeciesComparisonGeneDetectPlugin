#include "SettingsAction.h"
#include "CrossSpeciesComparisonGeneDetectPlugin.h"
#include<string>  
#include <QFileDialog>
#include <QPageLayout>
#include <QWebEngineView>
using namespace mv;
using namespace mv::gui;

SettingsAction::SettingsAction(CrossSpeciesComparisonGeneDetectPlugin& CrossSpeciesComparisonGeneDetectPlugin) :
    WidgetAction(&CrossSpeciesComparisonGeneDetectPlugin, "CrossSpeciesComparisonGeneDetectPlugin"),
    _crossSpeciesComparisonGeneDetectPlugin(CrossSpeciesComparisonGeneDetectPlugin),
    _tableModel(this, "Table Model"),
    _selectedGene(this, "Selected Gene"),
    _optionSelectionAction(*this)
{
    setSerializationName("CrossSpeciesComparisonGeneDetectPlugin Settings");
    _tableModel.setSerializationName("CSCGDV::Table Model");
    _selectedGene.setSerializationName("CSCGDV::Selected Gene");

    _selectedGene.setDisabled(true);
    _selectedGene.setString("");

}


SettingsAction::Widget::Widget(QWidget* parent, SettingsAction* SettingsAction) :
    WidgetActionWidget(parent, SettingsAction)
{
}



SettingsAction::OptionSelectionAction::Widget::Widget(QWidget* parent, OptionSelectionAction* optionSelectionAction) :
    WidgetActionWidget(parent, optionSelectionAction)
{
    //auto& SettingsAction = optionSelectionAction->_settingsAction;

    //
    //auto selectiondeStats1Widget = SettingsAction._deStatsDataset1Action.createWidget(this);
    //selectiondeStats1Widget->findChild<QComboBox*>("ComboBox")->setSizeAdjustPolicy(QComboBox::AdjustToContents);
 

    //auto selectiondeStats2Widget = SettingsAction._deStatsDataset2Action.createWidget(this);
    //selectiondeStats2Widget->findChild<QComboBox*>("ComboBox")->setSizeAdjustPolicy(QComboBox::AdjustToContents);

    //auto selectionExampledeStatsOptionLayout = new QFormLayout();
    //selectionExampledeStatsOptionLayout->setContentsMargins(0, 0, 0, 0);

    //selectionExampledeStatsOptionLayout->addRow(SettingsAction._deStatsDataset1Action.createLabelWidget(this), selectiondeStats1Widget);

    //selectionExampledeStatsOptionLayout->addRow(SettingsAction._deStatsDataset2Action.createLabelWidget(this), selectiondeStats2Widget);

    //setLayout(selectionExampledeStatsOptionLayout);
}

inline SettingsAction::OptionSelectionAction::OptionSelectionAction(SettingsAction& SettingsAction) :
    GroupAction(nullptr, "CrossSpeciesComparisonGeneDetectPluginOptionSelectionAction"),
    _settingsAction(SettingsAction)
{
    setText("Options");
    setIcon(Application::getIconFont("FontAwesome").getIcon("database"));
    addAction(&_settingsAction.getTableModelAction());
    addAction(&_settingsAction.getSelectedGeneAction());
}


void SettingsAction::fromVariantMap(const QVariantMap& variantMap)
{
    WidgetAction::fromVariantMap(variantMap);

    _tableModel.fromParentVariantMap(variantMap);
    _selectedGene.fromParentVariantMap(variantMap);
}

QVariantMap SettingsAction::toVariantMap() const
{
    QVariantMap variantMap = WidgetAction::toVariantMap();

    _tableModel.insertIntoVariantMap(variantMap);
    _selectedGene.insertIntoVariantMap(variantMap);

    return variantMap;
}
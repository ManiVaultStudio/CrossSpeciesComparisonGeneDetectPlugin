#include "SettingsAction.h"
#include "CrossSpeciesComparisonGeneDetectPlugin.h"
#include<string>  
#include <QFileDialog>
#include <QPageLayout>
#include <QWebEngineView>
using namespace mv;
using namespace mv::gui;

SettingsAction::SettingsAction(CrossSpeciesComparisonGeneDetectPlugin& CrossSpeciesComparisonGeneDetectPlugin) :
    WidgetAction(&CrossSpeciesComparisonGeneDetectPlugin, "CrossSpeciesComparisonGeneDetectPlugin Settings"),
    _crossSpeciesComparisonGeneDetectPlugin(CrossSpeciesComparisonGeneDetectPlugin),
    _tableModel(this, "Table Model"),
    _selectedGene(this, "Selected Gene"),
    _treeDataset(this, "Tree Dataset"),
    _selectedRowIndex(this, "Selected Row Index"),
    _optionSelectionAction(*this)
{
    setSerializationName("CSCGDV:CrossSpeciesComparison Gene Detect Plugin Settings");
    _tableModel.setSerializationName("CSCGDV:Table Model");
    _selectedGene.setSerializationName("CSCGDV:Selected Gene");

    _selectedGene.setDisabled(true);
    _selectedGene.setString("");

    _treeDataset.setSerializationName("CSCGDV:Tree Dataset");
    _selectedRowIndex.setSerializationName("CSCGDV:Selected Row Index");
    _selectedRowIndex.setDisabled(true);
    _selectedRowIndex.setString("");

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
    addAction(&_settingsAction.getTableModelAction());
    addAction(&_settingsAction.getSelectedGeneAction());
    addAction(&_settingsAction.getTreeDatasetAction());
}


void SettingsAction::fromVariantMap(const QVariantMap& variantMap)
{
    WidgetAction::fromVariantMap(variantMap);

    _tableModel.fromParentVariantMap(variantMap);
    _selectedGene.fromParentVariantMap(variantMap);
    _treeDataset.fromParentVariantMap(variantMap);
    _selectedRowIndex.fromParentVariantMap(variantMap);
}

QVariantMap SettingsAction::toVariantMap() const
{
    QVariantMap variantMap = WidgetAction::toVariantMap();

    _tableModel.insertIntoVariantMap(variantMap);
    _selectedGene.insertIntoVariantMap(variantMap);
    _treeDataset.insertIntoVariantMap(variantMap);
    _selectedRowIndex.insertIntoVariantMap(variantMap);

    return variantMap;
}
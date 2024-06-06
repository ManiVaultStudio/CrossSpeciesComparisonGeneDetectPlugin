#pragma once
#include <actions/WidgetAction.h>
#include <actions/IntegralAction.h>
#include <actions/OptionAction.h>
#include <actions/OptionsAction.h>
#include <actions/ToggleAction.h>
#include "actions/DatasetPickerAction.h"
#include "PointData/PointData.h"
#include "ClusterData/ClusterData.h"
#include "event/EventListener.h"
#include "actions/Actions.h"
#include "Plugin.h"
#include "DataHierarchyItem.h"
#include "Set.h"
#include <AnalysisPlugin.h>
#include <memory>
#include <algorithm>    
#include <QDebug>
#include <QLabel>
#include <QComboBox>
#include <QGroupBox>
#include <QPushButton>
#include <QGridLayout>
#include <QFormLayout>
#include <QString>
#include <string>
#include <event/Event.h>
#include <QDebug>
#include <QLabel>
#include <string>
#include "actions/VariantAction.h"
#include "actions/GroupAction.h"
using namespace mv::gui;
class QMenu;
class CrossSpeciesComparisonGeneDetectPlugin;

class FetchMetaData;
namespace mv
{
    class CoreInterface;
}

class SettingsAction : public WidgetAction
{
public:
    class OptionSelectionAction : public GroupAction
    {
    protected:
        class Widget : public mv::gui::WidgetActionWidget {
        public:
            Widget(QWidget* parent, OptionSelectionAction* optionSelectionAction);
        };

        QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
            return new OptionSelectionAction::Widget(parent, this);
        };

    public:
        OptionSelectionAction(SettingsAction& SettingsAction);

    protected:
        SettingsAction& _settingsAction;
    };



protected:

    class Widget : public mv::gui::WidgetActionWidget {
    public:
        Widget(QWidget* parent, SettingsAction* SettingsAction);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new SettingsAction::Widget(parent, this);
    };

public:
    SettingsAction(CrossSpeciesComparisonGeneDetectPlugin& CrossSpeciesComparisonGeneDetectPlugins);

public: // Action getters

    VariantAction& getTableModelAction() { return _tableModel; }
    StringAction& getSelectedGeneAction() { return _selectedGene; }
    StringAction&  getSelectedRowIndexAction() { return _selectedRowIndex; }
    DatasetPickerAction& getTreeDatasetAction() { return _treeDataset; }
    OptionSelectionAction& getOptionSelectionAction() { return _optionSelectionAction; }


public: // Serialization

    /**
     * Load widget action from variant map
     * @param Variant map representation of the widget action
     */
    void fromVariantMap(const QVariantMap& variantMap) override;

    /**
     * Save widget action to variant map
     * @return Variant map representation of the widget action
     */
    QVariantMap toVariantMap() const override;

protected:
    CrossSpeciesComparisonGeneDetectPlugin& _crossSpeciesComparisonGeneDetectPlugin;
    VariantAction                 _tableModel;
    StringAction                  _selectedGene;
    DatasetPickerAction          _treeDataset;
    StringAction                _selectedRowIndex;
    OptionSelectionAction         _optionSelectionAction;
    
};
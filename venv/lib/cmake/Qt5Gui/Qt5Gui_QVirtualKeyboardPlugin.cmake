
add_library(Qt5::QVirtualKeyboardPlugin MODULE IMPORTED)

_populate_Gui_plugin_properties(QVirtualKeyboardPlugin RELEASE "platforminputcontexts/libqtvirtualkeyboardplugin.dylib")

list(APPEND Qt5Gui_PLUGINS Qt5::QVirtualKeyboardPlugin)

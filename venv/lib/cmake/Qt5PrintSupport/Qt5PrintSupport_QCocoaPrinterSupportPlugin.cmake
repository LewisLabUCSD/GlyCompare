
add_library(Qt5::QCocoaPrinterSupportPlugin MODULE IMPORTED)

_populate_PrintSupport_plugin_properties(QCocoaPrinterSupportPlugin RELEASE "printsupport/libcocoaprintersupport.dylib")

list(APPEND Qt5PrintSupport_PLUGINS Qt5::QCocoaPrinterSupportPlugin)

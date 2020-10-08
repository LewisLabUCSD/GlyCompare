
add_library(Qt5::GeoServiceProviderFactoryEsri MODULE IMPORTED)

_populate_Location_plugin_properties(GeoServiceProviderFactoryEsri RELEASE "geoservices/libqtgeoservices_esri.dylib")

list(APPEND Qt5Location_PLUGINS Qt5::GeoServiceProviderFactoryEsri)

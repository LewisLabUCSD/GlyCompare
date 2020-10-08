
add_library(Qt5::CoreAudioPlugin MODULE IMPORTED)

_populate_Multimedia_plugin_properties(CoreAudioPlugin RELEASE "audio/libqtaudio_coreaudio.dylib")

list(APPEND Qt5Multimedia_PLUGINS Qt5::CoreAudioPlugin)

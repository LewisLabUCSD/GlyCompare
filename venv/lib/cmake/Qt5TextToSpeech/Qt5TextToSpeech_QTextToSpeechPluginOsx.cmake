
add_library(Qt5::QTextToSpeechPluginOsx MODULE IMPORTED)

_populate_TextToSpeech_plugin_properties(QTextToSpeechPluginOsx RELEASE "texttospeech/libqtexttospeech_speechosx.dylib")

list(APPEND Qt5TextToSpeech_PLUGINS Qt5::QTextToSpeechPluginOsx)

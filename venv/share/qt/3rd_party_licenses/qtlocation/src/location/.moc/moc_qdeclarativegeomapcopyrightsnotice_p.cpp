/****************************************************************************
** Meta object code from reading C++ file 'qdeclarativegeomapcopyrightsnotice_p.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.12.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../declarativemaps/qdeclarativegeomapcopyrightsnotice_p.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'qdeclarativegeomapcopyrightsnotice_p.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.12.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_QDeclarativeGeoMapCopyrightNotice_t {
    QByteArrayData data[16];
    char stringdata0[263];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_QDeclarativeGeoMapCopyrightNotice_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_QDeclarativeGeoMapCopyrightNotice_t qt_meta_stringdata_QDeclarativeGeoMapCopyrightNotice = {
    {
QT_MOC_LITERAL(0, 0, 33), // "QDeclarativeGeoMapCopyrightNo..."
QT_MOC_LITERAL(1, 34, 13), // "linkActivated"
QT_MOC_LITERAL(2, 48, 0), // ""
QT_MOC_LITERAL(3, 49, 4), // "link"
QT_MOC_LITERAL(4, 54, 16), // "mapSourceChanged"
QT_MOC_LITERAL(5, 71, 22), // "backgroundColorChanged"
QT_MOC_LITERAL(6, 94, 5), // "color"
QT_MOC_LITERAL(7, 100, 17), // "styleSheetChanged"
QT_MOC_LITERAL(8, 118, 10), // "styleSheet"
QT_MOC_LITERAL(9, 129, 24), // "copyrightsVisibleChanged"
QT_MOC_LITERAL(10, 154, 17), // "copyrightsChanged"
QT_MOC_LITERAL(11, 172, 15), // "copyrightsImage"
QT_MOC_LITERAL(12, 188, 14), // "copyrightsHtml"
QT_MOC_LITERAL(13, 203, 29), // "onCopyrightsStyleSheetChanged"
QT_MOC_LITERAL(14, 233, 9), // "mapSource"
QT_MOC_LITERAL(15, 243, 19) // "QDeclarativeGeoMap*"

    },
    "QDeclarativeGeoMapCopyrightNotice\0"
    "linkActivated\0\0link\0mapSourceChanged\0"
    "backgroundColorChanged\0color\0"
    "styleSheetChanged\0styleSheet\0"
    "copyrightsVisibleChanged\0copyrightsChanged\0"
    "copyrightsImage\0copyrightsHtml\0"
    "onCopyrightsStyleSheetChanged\0mapSource\0"
    "QDeclarativeGeoMap*"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_QDeclarativeGeoMapCopyrightNotice[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       2,   74, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       5,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   54,    2, 0x06 /* Public */,
       4,    0,   57,    2, 0x06 /* Public */,
       5,    1,   58,    2, 0x06 /* Public */,
       7,    1,   61,    2, 0x06 /* Public */,
       9,    0,   64,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      10,    1,   65,    2, 0x0a /* Public */,
      10,    1,   68,    2, 0x0a /* Public */,
      13,    1,   71,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QColor,    6,
    QMetaType::Void, QMetaType::QString,    8,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void, QMetaType::QImage,   11,
    QMetaType::Void, QMetaType::QString,   12,
    QMetaType::Void, QMetaType::QString,    8,

 // properties: name, type, flags
      14, 0x80000000 | 15, 0x0049510b,
       8, QMetaType::QString, 0x00495103,

 // properties: notify_signal_id
       1,
       3,

       0        // eod
};

void QDeclarativeGeoMapCopyrightNotice::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<QDeclarativeGeoMapCopyrightNotice *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->linkActivated((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: _t->mapSourceChanged(); break;
        case 2: _t->backgroundColorChanged((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 3: _t->styleSheetChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: _t->copyrightsVisibleChanged(); break;
        case 5: _t->copyrightsChanged((*reinterpret_cast< const QImage(*)>(_a[1]))); break;
        case 6: _t->copyrightsChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 7: _t->onCopyrightsStyleSheetChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (QDeclarativeGeoMapCopyrightNotice::*)(const QString & );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&QDeclarativeGeoMapCopyrightNotice::linkActivated)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (QDeclarativeGeoMapCopyrightNotice::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&QDeclarativeGeoMapCopyrightNotice::mapSourceChanged)) {
                *result = 1;
                return;
            }
        }
        {
            using _t = void (QDeclarativeGeoMapCopyrightNotice::*)(const QColor & );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&QDeclarativeGeoMapCopyrightNotice::backgroundColorChanged)) {
                *result = 2;
                return;
            }
        }
        {
            using _t = void (QDeclarativeGeoMapCopyrightNotice::*)(const QString & );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&QDeclarativeGeoMapCopyrightNotice::styleSheetChanged)) {
                *result = 3;
                return;
            }
        }
        {
            using _t = void (QDeclarativeGeoMapCopyrightNotice::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&QDeclarativeGeoMapCopyrightNotice::copyrightsVisibleChanged)) {
                *result = 4;
                return;
            }
        }
    }
#ifndef QT_NO_PROPERTIES
    else if (_c == QMetaObject::ReadProperty) {
        auto *_t = static_cast<QDeclarativeGeoMapCopyrightNotice *>(_o);
        Q_UNUSED(_t)
        void *_v = _a[0];
        switch (_id) {
        case 0: *reinterpret_cast< QDeclarativeGeoMap**>(_v) = _t->mapSource(); break;
        case 1: *reinterpret_cast< QString*>(_v) = _t->styleSheet(); break;
        default: break;
        }
    } else if (_c == QMetaObject::WriteProperty) {
        auto *_t = static_cast<QDeclarativeGeoMapCopyrightNotice *>(_o);
        Q_UNUSED(_t)
        void *_v = _a[0];
        switch (_id) {
        case 0: _t->setMapSource(*reinterpret_cast< QDeclarativeGeoMap**>(_v)); break;
        case 1: _t->setStyleSheet(*reinterpret_cast< QString*>(_v)); break;
        default: break;
        }
    } else if (_c == QMetaObject::ResetProperty) {
    }
#endif // QT_NO_PROPERTIES
}

QT_INIT_METAOBJECT const QMetaObject QDeclarativeGeoMapCopyrightNotice::staticMetaObject = { {
    &QQuickPaintedItem::staticMetaObject,
    qt_meta_stringdata_QDeclarativeGeoMapCopyrightNotice.data,
    qt_meta_data_QDeclarativeGeoMapCopyrightNotice,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *QDeclarativeGeoMapCopyrightNotice::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *QDeclarativeGeoMapCopyrightNotice::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_QDeclarativeGeoMapCopyrightNotice.stringdata0))
        return static_cast<void*>(this);
    return QQuickPaintedItem::qt_metacast(_clname);
}

int QDeclarativeGeoMapCopyrightNotice::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QQuickPaintedItem::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 8)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 8;
    }
#ifndef QT_NO_PROPERTIES
    else if (_c == QMetaObject::ReadProperty || _c == QMetaObject::WriteProperty
            || _c == QMetaObject::ResetProperty || _c == QMetaObject::RegisterPropertyMetaType) {
        qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    } else if (_c == QMetaObject::QueryPropertyDesignable) {
        _id -= 2;
    } else if (_c == QMetaObject::QueryPropertyScriptable) {
        _id -= 2;
    } else if (_c == QMetaObject::QueryPropertyStored) {
        _id -= 2;
    } else if (_c == QMetaObject::QueryPropertyEditable) {
        _id -= 2;
    } else if (_c == QMetaObject::QueryPropertyUser) {
        _id -= 2;
    }
#endif // QT_NO_PROPERTIES
    return _id;
}

// SIGNAL 0
void QDeclarativeGeoMapCopyrightNotice::linkActivated(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void QDeclarativeGeoMapCopyrightNotice::mapSourceChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}

// SIGNAL 2
void QDeclarativeGeoMapCopyrightNotice::backgroundColorChanged(const QColor & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void QDeclarativeGeoMapCopyrightNotice::styleSheetChanged(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void QDeclarativeGeoMapCopyrightNotice::copyrightsVisibleChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 4, nullptr);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE

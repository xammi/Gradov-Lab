#-------------------------------------------------
#
# Project created by QtCreator 2016-09-21T00:35:15
#
#-------------------------------------------------

QT       += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = LR_Methods
TEMPLATE = app
CONFIG += c++11


SOURCES += main.cpp\
    interpolation.cpp \
    gauss.cpp \
    mainwindow.cpp

HEADERS  += mainwindow.h\
            interpolation.h \
    resources.h \
    gauss.h

FORMS    += mainwindow.ui

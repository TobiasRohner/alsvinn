#include "alsfvm/init/PythonInitialData.hpp"
#ifdef _DEBUG
#undef _DEBUG
#include <Python.h>
#define _DEBUG
#else
#include <Python.h>
#endif
#include <iostream>
#include <sstream>
#include "alsfvm/volume/volume_foreach.hpp"
#include "alsfvm/python/PythonInterpreter.hpp"
#define L std::cout << __FILE__ << ":" << __LINE__ << std::endl;
using namespace alsfvm::python;
namespace alsfvm { namespace init {

namespace {


    // Adds indent to each line, eg transfer
    // "a\nb\n" to
    // "    a\n    b\n"
    void addIndent(const std::string& inputString, std::stringstream& output) {
        std::istringstream userSpecifiedFunction(inputString);

        std::string line;
        while (std::getline(userSpecifiedFunction, line)) {
            output <<"    "<< line << std::endl;
        }
    }
}
PythonInitialData::PythonInitialData(const std::string &programString)
    : programString(programString)
{

}

void PythonInitialData::setInitialData(volume::Volume& conservedVolume,
                                       volume::Volume& extraVolume,
                                       volume::Volume& primitiveVolume,
                                       equation::CellComputer& cellComputer,
                                       grid::Grid& grid)
{

    auto& interpreter = PythonInterpreter::getInstance();

    // We need to add namespaces for the function to live in, this way
    // we can actually reference later.


    auto globalNamespace = interpreter.getGlobalNamespace();
    PyObject* moduleLocal(PyImport_AddModule("alsfvm"));
    PythonObjectHolder localNamespace(PyModule_GetDict(moduleLocal));
    if (PyErr_Occurred()) {
        PyErr_Print();
        THROW("Error in python script.");
    }
    // This will hold the inputs for our function. We allocate this once,
    // and use it several times.
    PythonObjectHolder argumentTuple(PyTuple_New(4));


    // Now we declare the wrappers around the function.
    std::stringstream functionStringStream;

    functionStringStream << "def initial_data(x, y, z, output):\n";
    addIndent("from math import *", functionStringStream);

    // Now we need to add the variables we need to write:
    for(size_t i = 0; i < primitiveVolume.getNumberOfVariables(); ++i) {
        // We set them to None, that way we can check at the end if they are checked.
        addIndent(primitiveVolume.getName(i) + " = 0.0", functionStringStream);
    }

    addIndent(programString, functionStringStream);

    // Add code to store the variables:
    for(size_t i = 0; i < primitiveVolume.getNumberOfVariables(); ++i) {
        // We set them to None, that way we can check at the end if they are checked.
        const auto& name = primitiveVolume.getName(i);
        addIndent(std::string("output['") + name + "'] = " + name, functionStringStream);
    }


    PyRun_String(functionStringStream.str().c_str(),
                 Py_file_input, globalNamespace, localNamespace.object);


    if (PyErr_Occurred()) {
        PyErr_Print();
        THROW("Error in python script.");
    }


    PythonObjectHolder initialValueFunction(PyObject_GetAttrString(moduleLocal, "initial_data"));
    if ( PyErr_Occurred()) {
        PyErr_Print();
        THROW("Python error occured");
    }

    // loop through the map and set the initial values
    volume::for_each_midpoint(primitiveVolume, grid,
                              [&](real x, real y, real z, size_t index) {

		if (PyErr_Occurred()) {
			PyErr_Print();
			THROW("Python error occured");
		}
        PyObject* outputMap(PyDict_New());

        PyObject* xObject(PyFloat_FromDouble(x));
        PyObject* yObject(PyFloat_FromDouble(y));
        PyObject* zObject(PyFloat_FromDouble(z));
        if ( PyErr_Occurred()) {
            PyErr_Print();
            THROW("Python error occured");
        }
        PyTuple_SetItem(argumentTuple.object, 0, xObject);
        PyTuple_SetItem(argumentTuple.object, 1, yObject);
        PyTuple_SetItem(argumentTuple.object, 2, zObject);
        PyTuple_SetItem(argumentTuple.object, 3, outputMap);

        if ( PyErr_Occurred()) {
            PyErr_Print();
            THROW("Python error occured");
        }
        PyObject_CallObject(initialValueFunction.object,
                                               argumentTuple.object);
		if (PyErr_Occurred()) {
			PyErr_Print();
			THROW("Python error occured");
		}
        // Loop through each variable and set it in the primitive variables:
        for(size_t var = 0; var <  primitiveVolume.getNumberOfVariables(); ++var) {
            const auto& name = primitiveVolume.getName(var);
            PyObject* floatObject(PyDict_GetItemString(outputMap, name.c_str()));
            const double value = PyFloat_AsDouble(floatObject);
            primitiveVolume.getScalarMemoryArea(var)->getPointer()[index] = value;
        }
		if (PyErr_Occurred()) {
			PyErr_Print();
			THROW("Python error occured");
		}
    });


    if (PyErr_Occurred()) {
        PyErr_Print();
        THROW("Error in python script.");
    }

    cellComputer.computeFromPrimitive(primitiveVolume, conservedVolume, extraVolume);
}



}
}

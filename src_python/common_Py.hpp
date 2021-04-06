#ifndef common_Py_h
#define common_Py_h

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <memory>
#include <array>

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include<string>
#include "ndarrayobject.h"
#include "float.h"
#include "console.hpp"
#include "QCM.hpp"
#include "qcm_ED.hpp"
#include "parameter_set.hpp"
#include "lattice_model.hpp"

extern shared_ptr<parameter_set> param_set;
extern shared_ptr<lattice_model> qcm_model;
extern vector<string> target_sectors;

void check_signals();
string Py2string(PyObject* PyObj);
vector3D<int64_t> position_from_Py(PyArrayObject *k_pyobj);
vector3D<double> wavevector_from_Py(PyArrayObject *k_pyobj);
vector<vector3D<double>> wavevectors_from_Py(PyArrayObject *k_pyobj);
vector<vector3D<int64_t>> intvec3D_from_Py(PyArrayObject *k_pyobj);
vector<int64_t> intvector_from_Py(PyArrayObject *k_pyobj);
vector<double> doublematrix_from_Py(PyArrayObject *k_pyobj);
vector<string> strings_from_PyList(PyObject* lst);
vector<double> doubles_from_Py(PyObject* lst);

map<string,double> py_dict_to_map(PyObject *D);
vector<vector<int>> intmatrix_from_Py(PyArrayObject *k_pyobj);
vector<vector<double>> pos_from_Py(PyArrayObject *k_pyobj);

#endif
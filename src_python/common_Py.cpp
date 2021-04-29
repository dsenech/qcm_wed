/**
 This files defines the Python wrappers and the python module
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL QCM_ARRAY_API
#define NO_IMPORT_ARRAY

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

//==============================================================================
// python utilities

void check_signals()
{
  if(PyErr_CheckSignals() == -1) throw;
}
//------------------------------------------------------------------------------
string Py2string(PyObject* PyObj)
{
  Py_ssize_t s = PyUnicode_GetLength(PyObj);
  const char* c = PyUnicode_AsUTF8(PyObj);
  return string(c, s);
}
//------------------------------------------------------------------------------
vector3D<int64_t> position_from_Py(PyArrayObject *k_pyobj)
{
  vector3D<int64_t> k;
  if(PyArray_Check(k_pyobj)){
    auto dims = PyArray_DIMS(k_pyobj);
    size_t nc = dims[0];
    if(nc!=3) qcm_throw("The supplied wavevector must have 3 components!");
    k.x = *(int*)PyArray_GETPTR1(k_pyobj, 0);
    k.y = *(int*)PyArray_GETPTR1(k_pyobj, 1);
    k.z = *(int*)PyArray_GETPTR1(k_pyobj, 2);
  }
  else if(PyList_Check(k_pyobj)){
    PyObject* PyObj = (PyObject*)k_pyobj;
    if(PyList_Size(PyObj) != 3) qcm_throw("a position must have 3 components");
    PyObject* pkey;
    pkey = PyList_GetItem(PyObj,0);
    if(!PyLong_Check(pkey)) qcm_throw("element no 1 of list is not an integer");
    k.x = PyLong_AsLong(pkey);
    pkey = PyList_GetItem(PyObj,1);
    if(!PyLong_Check(pkey)) qcm_throw("element no 2 of list is not an integer");
    k.y = PyLong_AsLong(pkey);
    pkey = PyList_GetItem(PyObj,2);
    if(!PyLong_Check(pkey)) qcm_throw("element no 3 of list is not an integer");
    k.z = PyLong_AsLong(pkey);
  }
  else qcm_throw("object is neither a numpy array or a python list");
  return k;
}
//------------------------------------------------------------------------------
vector3D<double> wavevector_from_Py(PyArrayObject *k_pyobj)
{
  vector3D<double> k;
  if(PyArray_Check(k_pyobj)){
    auto dims = PyArray_DIMS(k_pyobj);
    size_t nc = dims[0];
    if(nc!=3) qcm_throw("The supplied wavevector must have 3 components!");
    k.x = *(double*)PyArray_GETPTR1(k_pyobj, 0);
    k.y = *(double*)PyArray_GETPTR1(k_pyobj, 1);
    k.z = *(double*)PyArray_GETPTR1(k_pyobj, 2);
  }
  else if(PyList_Check(k_pyobj)){
    PyObject* PyObj = (PyObject*)k_pyobj;
    if(PyList_Size(PyObj) != 3) qcm_throw("a wavevector must have 3 components");
    k.x = PyFloat_AsDouble(PyList_GetItem(PyObj,0));
    k.y = PyFloat_AsDouble(PyList_GetItem(PyObj,1));
    k.z = PyFloat_AsDouble(PyList_GetItem(PyObj,2));
  }
  else qcm_throw("object is neither a numpy array or a python list");
  return k;
}
//------------------------------------------------------------------------------
vector<vector3D<double>> wavevectors_from_Py(PyArrayObject *k_pyobj)
{
  if(PyArray_NDIM(k_pyobj) != 2) qcm_throw("The supplied wavevectors must form a 2D Numpy array!");
  
  vector<npy_intp> dims(2);
  dims[0] = *PyArray_DIMS(k_pyobj);
  dims[1] = *(PyArray_DIMS(k_pyobj)+1);
  
  if(dims[1] != 3) qcm_throw("The supplied wavevectors must have 3 components!");
  vector<vector3D<double>> k(dims[0]);
  
  for(size_t i=0; i<k.size(); i++){
    k[i].x = *(double*)PyArray_GETPTR2(k_pyobj, i, 0);
    k[i].y = *(double*)PyArray_GETPTR2(k_pyobj, i, 1);
    k[i].z = *(double*)PyArray_GETPTR2(k_pyobj, i, 2);
  }
  return k;
}
//------------------------------------------------------------------------------
vector<vector3D<int64_t>> intvec3D_from_Py(PyArrayObject *k_pyobj)
{
  vector<vector3D<int64_t>> k;
  if(PyArray_Check(k_pyobj)){
    if(PyArray_NDIM(k_pyobj) != 2) qcm_throw("The supplied wavevectors must form a 2D Numpy array!");
    
    vector<npy_intp> dims(2);
    dims[0] = *PyArray_DIMS(k_pyobj);
    dims[1] = *(PyArray_DIMS(k_pyobj)+1);
    
    if(dims[1] != 3) qcm_throw("The supplied integer vectors must have 3 components!");
    k.assign(dims[0], vector3D<int64_t>());
    
    for(size_t i=0; i<k.size(); i++){
      k[i].x = *(int64_t*)PyArray_GETPTR2(k_pyobj, i, 0);
      k[i].y = *(int64_t*)PyArray_GETPTR2(k_pyobj, i, 1);
      k[i].z = *(int64_t*)PyArray_GETPTR2(k_pyobj, i, 2);
    }
  }
  else if(PyList_Check(k_pyobj)){
    PyObject* PyObj = (PyObject*)k_pyobj;
    size_t n = PyList_Size(PyObj);
    k.assign(n, vector3D<int64_t>());
    for(size_t i=0; i<k.size(); i++){
      PyObject* PyObj2 = PyList_GetItem(PyObj,i);
      if(PyList_Size(PyObj2) != 3) qcm_throw("a position must have 3 components");
      k[i].x = PyLong_AsLong(PyList_GetItem(PyObj2,0));
      k[i].y = PyLong_AsLong(PyList_GetItem(PyObj2,1));
      k[i].z = PyLong_AsLong(PyList_GetItem(PyObj2,2));
    }
  }
  return k;
}
//------------------------------------------------------------------------------
vector<int64_t> intvector_from_Py(PyArrayObject *k_pyobj)
{
  vector<int64_t> k;
  if(PyArray_Check(k_pyobj)){
    vector<npy_intp> dims(2);
    dims[0] = *PyArray_DIMS(k_pyobj);
    dims[1] = *(PyArray_DIMS(k_pyobj)+1);
    
    if(dims[1] != 3) qcm_throw("The supplied integer vectors must have 3 components!");
    k.assign(3*dims[0], 0);
    
    int l=0;
    for(size_t i=0; i<dims[0]; i++){
      for(size_t j=0; j<3; j++){
        k[l++] = *(int64_t*)PyArray_GETPTR2(k_pyobj, i, j);
      }
    }
  }
  else if(PyList_Check(k_pyobj)){
    PyObject* PyObj = (PyObject*)k_pyobj;
    size_t n = PyList_Size(PyObj);
    k.assign(3*n, 0);
    int l=0;
    for(size_t i=0; i<n; i++){
      PyObject* PyObj2 = PyList_GetItem(PyObj,i);
      if(PyList_Size(PyObj2) != 3) qcm_throw("the integer vectors must have 3 components");
      k[l++] = PyLong_AsLong(PyList_GetItem(PyObj2,0));
      k[l++] = PyLong_AsLong(PyList_GetItem(PyObj2,1));
      k[l++] = PyLong_AsLong(PyList_GetItem(PyObj2,2));
    }
  }
  return k;
}
//------------------------------------------------------------------------------
vector<double> doublematrix_from_Py(PyArrayObject *k_pyobj)
{
  vector<double> k;

  if(PyArray_Check(k_pyobj)){
    if(PyArray_NDIM(k_pyobj) != 2) qcm_throw("The supplied wavevectors must form a 2D Numpy array!");
    
    vector<npy_intp> dims(2);
    dims[0] = *PyArray_DIMS(k_pyobj);
    dims[1] = *(PyArray_DIMS(k_pyobj)+1);
    
    if(dims[1] != 3) qcm_throw("The supplied integer vectors must have 3 components!");
    k.assign(3*dims[0], 0.0);
    
    int l=0;
    for(size_t i=0; i<dims[0]; i++){
      for(size_t j=0; j<3; j++){
        k[l++] = *(int64_t*)PyArray_GETPTR2(k_pyobj, i, j);
      }
    }
  }
  else if(PyList_Check(k_pyobj)){
    PyObject* PyObj = (PyObject*)k_pyobj;
    size_t n = PyList_Size(PyObj);
    k.assign(3*n, 0);
    int l=0;
    for(size_t i=0; i<n; i++){
      PyObject* PyObj2 = PyList_GetItem(PyObj,i);
      if(PyList_Size(PyObj2) != 3) qcm_throw("the vectors must have 3 components");
      k[l++] = PyFloat_AsDouble(PyList_GetItem(PyObj2,0));
      k[l++] = PyFloat_AsDouble(PyList_GetItem(PyObj2,1));
      k[l++] = PyFloat_AsDouble(PyList_GetItem(PyObj2,2));
    }
  }
  return k;
}
//------------------------------------------------------------------------------
vector<string> strings_from_PyList(PyObject* lst)
{
  if(!PyList_Check(lst)) qcm_throw("expected a list of strings");
  size_t n = PyList_Size(lst);
  vector<string> out(n);
  for(size_t i=0; i<n; i++){
    PyObject* pkey = PyList_GetItem(lst,i);
    Py_ssize_t s = PyUnicode_GetLength(pkey);
    const char* key_char = PyUnicode_AsUTF8(pkey);
    out[i] = string(key_char, s);
  }
  return out;
}
//------------------------------------------------------------------------------
vector<double> doubles_from_Py(PyObject* lst)
{
  vector<double> out;
  if(PyList_Check(lst)){
    size_t n = PyList_Size(lst);
    out.resize(n);
    for(size_t i=0; i<n; i++){
      out[i] = PyFloat_AsDouble(PyList_GetItem(lst,i));
    }
  }
  else if(PyArray_Check(lst)){
    PyArrayObject* alst = (PyArrayObject*)(lst);
    if(PyArray_NDIM(alst) != 1) qcm_throw("The supplied NumPy array must be one-dimensional!");
    vector<npy_intp> dims(1);
    dims[0] = *PyArray_DIMS(alst);
    out.resize(dims[0]); 
    int l=0;
    for(size_t i=0; i<dims[0]; i++){
      out[i] = *(double*)PyArray_GETPTR1(alst, i);
    }
  }
  else{
    qcm_throw("expected list or NumPy array and got something else!");
  }
  return out;
}

//------------------------------------------------------------------------------
map<string,double> py_dict_to_map(PyObject *D){
  
  if(!PyDict_Check(D)) qcm_ED_throw("argument of dict_to_map() is not a python dictionary");

  map<string,double> the_map;
  Py_ssize_t ppos=0;
  PyObject *key = nullptr;
  PyObject *v = nullptr;
  while(PyDict_Next(D, &ppos, &key, &v)){
    Py_ssize_t s = PyUnicode_GetLength(key);
    const char* key_char = PyUnicode_AsUTF8(key);
    the_map[string(key_char, s)] = PyFloat_AsDouble(v);
  }
  return the_map;
}
//------------------------------------------------------------------------------
vector<vector<int>> intmatrix_from_Py(PyArrayObject *k_pyobj)
{
  vector<vector<int>> k;
  if(PyArray_Check(k_pyobj)){
    if(PyArray_NDIM(k_pyobj)!=2)
      qcm_ED_throw("generators or positions should form a two-dimensional array!");
    vector<npy_intp> dims(2);
    dims[0] = *PyArray_DIMS(k_pyobj);
    dims[1] = *(PyArray_DIMS(k_pyobj)+1);
    k.assign(dims[0], vector<int>(dims[1]));
    
    for(size_t i=0; i<dims[0]; i++){
      for(size_t j=0; j<dims[1]; j++){
        k[i][j]= *(int*)PyArray_GETPTR2(k_pyobj, i, j);
      }
    }
  }
  else if(PyList_Check(k_pyobj)){
    PyObject* PyObj = (PyObject*)k_pyobj;
    size_t n = PyList_Size(PyObj);
    k.assign(n, vector<int>());
    for(size_t i=0; i<n; i++){
      PyObject* PyObj2 = PyList_GetItem(PyObj,i);
      if(!PyList_Check(PyObj2))
        qcm_ED_throw("generator or position "+to_string(i+1)+" should be a list!");
      size_t m = PyList_Size(PyObj2);
      k[i].assign(m, 0);
      for(size_t j=0; j<m; j++){
        k[i][j] = (int)PyLong_AsLong(PyList_GetItem(PyObj2,j));
      }
    }
  }
  return k;
}

//------------------------------------------------------------------------------
vector<vector<double>> pos_from_Py(PyArrayObject *k_pyobj)
{
  vector<vector<double>> k;
  if(PyArray_Check(k_pyobj)){
    if(PyArray_NDIM(k_pyobj)!=2)
      qcm_ED_throw("positions should form a two-dimensional array!");
    vector<npy_intp> dims(2);
    dims[0] = *PyArray_DIMS(k_pyobj);
    dims[1] = *(PyArray_DIMS(k_pyobj)+1);
    
    k.assign(dims[0], vector<double>(dims[1]));
    
    for(size_t i=0; i<dims[0]; i++){
      for(size_t j=0; j<dims[1]; j++){
        k[i][j]= *(double*)PyArray_GETPTR2(k_pyobj, i, j);
      }
    }
  }
  else if(PyList_Check(k_pyobj)){
    PyObject* PyObj = (PyObject*)k_pyobj;
    size_t n = PyList_Size(PyObj);
    k.assign(n, vector<double>());
    for(size_t i=0; i<n; i++){
      PyObject* PyObj2 = PyList_GetItem(PyObj,i);
      if(!PyList_Check(PyObj2))
        qcm_ED_throw("position "+to_string(i+1)+" should be a list!");
      size_t m = PyList_Size(PyObj2);
      k[i].assign(m, 0);
      for(size_t j=0; j<m; j++){
        k[i][j] = (double)PyFloat_AsDouble(PyList_GetItem(PyObj2,j));
      }
    }
  }
  return k;
}

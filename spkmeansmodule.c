#include <Python.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "spkmeans.h"

/* C Api Declarations: */
static double **create2DArrayFromPyObject(PyObject *data, int n, int d);
static PyObject *create2DPyObject(double **data, int n, int d);
#define verifyNotNULL(var) if((var)==NULL) {printf("An Error Has Occured"); exit(-1);}

/*--------------------------------- c-api -----------------------------------------------------*/

static PyObject *sideGoals_Capi(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d ;
    double **data_points;
    char *goal;

    if (!PyArg_ParseTuple(args, "iisO", &n, &d,&goal, &points)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }

    data_points = create2DArrayFromPyObject(points, n, d);

    printByGoal(data_points,d,n,goal);

    free2DMalloc(data_points, n);

    return Py_True;
}

static PyObject *getT_Capi(PyObject *self, PyObject *args) {
    PyObject *pyPoints;
    PyObject *pyT;
    int k, n, d;
    double **T, **points;
    if (!PyArg_ParseTuple(args, "iiiO", &k, &n, &d, &pyPoints)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }

    points = create2DArrayFromPyObject(pyPoints, n, d);

    T=getT(points,d,n,&k);
    pyT = create2DPyObject(T,n,k);

    free2DMalloc(points, n);
    free_contiguous_mat(T);
    return Py_BuildValue("[Oi]", pyT, k);
}
static PyObject *kmeans_Capi(PyObject *self, PyObject *args) {
    PyObject *pyInitialCentroids, *pyPoints;
    int k, n, d, max_iter;
    double **initialCentroids, **points;
    if (!PyArg_ParseTuple(args, "iiiiOO", &k, &n, &d, &max_iter, &pyInitialCentroids, &pyPoints)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    initialCentroids = create2DArrayFromPyObject(pyInitialCentroids, k, d);
    points = create2DArrayFromPyObject(pyPoints, n, d);

    double **finalCentroids = kmeans(points, initialCentroids, k, d, n, max_iter);

    print_2d_array(finalCentroids, k, d);

    free2DMalloc(points, n);
    free2DMalloc(initialCentroids, k);

    return Py_True;
}
static double **create2DArrayFromPyObject(PyObject *data, int n, int d) {
    int i, j;
    double **points;
    PyObject *temp_point,*inner_item;

    points = (double **) malloc(n * sizeof(double *));
    verifyNotNULL(points)

    for (i = 0; i < n; i++) {
        double *vector = malloc(d * sizeof(double));
        verifyNotNULL(vector)


        temp_point = PyList_GetItem(data, i);
        for (j = 0; j < d; j++) {
            inner_item = PyList_GetItem(temp_point, j);
            vector[j] = PyFloat_AsDouble(inner_item);
        }
        points[i] = vector;
    }

    return points;
}
static PyObject *create2DPyObject(double** matrix, int n, int d) {
    int i,j;
    PyObject *currentVector, *pyMatrix,*num;
    pyMatrix = PyList_New(n);
    verifyNotNULL(pyMatrix)
    for (i = 0; i < n; i++) {
        currentVector = PyList_New(d);
        verifyNotNULL(currentVector)
        for (j = 0; j < d; j++) {
            num = PyFloat_FromDouble(matrix[i][j]);
            verifyNotNULL(num)
            PyList_SET_ITEM(currentVector, j, num);
        }
    PyList_SET_ITEM(pyMatrix, i, currentVector);
    }
    return pyMatrix;
}



static PyMethodDef capiMethods[] = {
        {"sideGoals",  (PyCFunction) sideGoals_Capi, METH_VARARGS, PyDoc_STR("Using C Extension to calculate and output a goal (wam,ddg,lnorm,jacobi).")},
        {"getT", (PyCFunction) getT_Capi, METH_VARARGS, PyDoc_STR("Calculate T.")},
        {"kmeans", (PyCFunction) kmeans_Capi, METH_VARARGS, PyDoc_STR("kmeans")},
        {NULL, NULL,                           0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,"myspkmeanssp",NULL,-1,capiMethods
};

PyMODINIT_FUNC
PyInit_myspkmeanssp(void)
{
    PyObject *n;
    n = PyModule_Create(&moduledef);
    if (!n)
        return NULL;
    return n;
}
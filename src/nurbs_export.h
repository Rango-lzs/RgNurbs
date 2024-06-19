#ifndef NURBS_EXPORT_HH
#define NURBS_EXPORT_HH

#define DLL_EXPORT __declspec(dllexport)
#define DLL_IMPORT __declspec(dllimport)

#ifdef RG_NURBS_MODULE
#define RG_API DLL_EXPORT
#else
#define RG_API DLL_IMPORT
#endif

#endif

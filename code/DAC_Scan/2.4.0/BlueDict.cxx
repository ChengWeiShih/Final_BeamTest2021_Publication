// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME BlueDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "Blue.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void delete_Blue(void *p);
   static void deleteArray_Blue(void *p);
   static void destruct_Blue(void *p);
   static void streamer_Blue(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Blue*)
   {
      ::Blue *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Blue >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Blue", ::Blue::Class_Version(), "Blue.h", 26,
                  typeid(::Blue), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Blue::Dictionary, isa_proxy, 16,
                  sizeof(::Blue) );
      instance.SetDelete(&delete_Blue);
      instance.SetDeleteArray(&deleteArray_Blue);
      instance.SetDestructor(&destruct_Blue);
      instance.SetStreamerFunc(&streamer_Blue);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Blue*)
   {
      return GenerateInitInstanceLocal((::Blue*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Blue*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Blue::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Blue::Class_Name()
{
   return "Blue";
}

//______________________________________________________________________________
const char *Blue::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Blue*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Blue::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Blue*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Blue::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Blue*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Blue::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Blue*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Blue::Streamer(TBuffer &R__b)
{
   // Stream an object of class Blue.

   TObject::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Blue(void *p) {
      delete ((::Blue*)p);
   }
   static void deleteArray_Blue(void *p) {
      delete [] ((::Blue*)p);
   }
   static void destruct_Blue(void *p) {
      typedef ::Blue current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Blue(TBuffer &buf, void *obj) {
      ((::Blue*)obj)->::Blue::Streamer(buf);
   }
} // end of namespace ROOT for class ::Blue

namespace {
  void TriggerDictionaryInitialization_BlueDict_Impl() {
    static const char* headers[] = {
"Blue.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/cwshih/.conda/envs/root_env/include/",
"/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/code/DAC_Scan/2.4.0/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "BlueDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$Blue.h")))  Blue;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "BlueDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "Blue.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Blue", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("BlueDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_BlueDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_BlueDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_BlueDict() {
  TriggerDictionaryInitialization_BlueDict_Impl();
}

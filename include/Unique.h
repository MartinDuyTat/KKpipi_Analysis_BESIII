// Copied from GGSZ code repository
#ifndef Unique_h
#define Unique_h

#include <iostream>
#include <string>
#include <unordered_set>

class Unique{
    /// Class that ensures RooFit objects have unique names
    /** If RooFit objects share a name, hell breaks loose
    due a lot of things being handled internally by name lookup.

    This class wraps around *anything* that can be constructed with
    a C string as the first argument, the name, and ensure nothing
    else has been created via GGSZUnique with the same name.

    If this is the case, _0 is added to the name until it is unique,
    and colourful warnings are thrown. This makes sure the code works,
    but is also annoying enough that one goes about fixing the naming
    issue :)

    The use is

      auto v = GGSZUnique::create<T>( args to T constructor );
      // v is a T*
      v->Print(); // shows what you just made!

    An example that will throw warnings

      auto v1 = GGSZUnique::create<RooRealVar>("a", "", 1, 0, 2);
      // ok!
      auto v2 = GGSZUnique::create<RooConstVar>("a", "", 4);
      // Warning! v2 is now names 'a_0'
      auto v3 = GGSZUnique::create<RooFormulaVar>("a", "@0+@1", 
                                                    RooArgList(*v1, *v2));
      // Warnings! v3 is now named 'a_0_0'

    */



public:

    // The public create function to be called
    // calls private template depending on whether T is pointer
    template <typename T, typename... Args>
    static T create (
      std::string name, 
      Args&&... args
      ) {
        return create(type<T>(), name, std::forward<Args>(args)...);
    } ;

private:

    // Private constructor ensures nonone instantiates GGSZUnique
    Unique(){};

    // Struct template
    template<typename T>
    struct type{};

    // Function that actually ensures a name is unique
    static std::string _make_unique(
      std::string& name
      ){
        static std::unordered_set<std::string> _names;
        auto success = _names.insert(name);
        if (!success.second){
            std::cerr << "A RooFit object with name '" << name << "' already exists! padding '_0' and trying again (FIX naming)\n";
            name += "_0";
            return _make_unique(name);
        }
        return name;
    };

    // Private create function for non-pointers
    template <typename T, typename... Args>
    static T create (
      type<T>,
      std::string& name,
      Args&&... args
      ) {
        name = _make_unique(name);
        T r(name.c_str(), std::forward<Args>(args)...);
        return r;
    }

    // Private create function for pointers
    template <typename T, typename... Args>
    static T* create (
      type<T*>,
      std::string& name,
      Args&&... args
      ) {
        name = _make_unique(name);
        T* r = new T(name.c_str(), std::forward<Args>(args)...);
        return r;
    }

} ;

#endif // GGSZUnique_h

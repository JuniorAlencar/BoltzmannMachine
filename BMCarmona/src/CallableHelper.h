#ifndef CALLABLE_HELPER_H
#define CALLABLE_HELPER_H

#include <utility>
#include <iostream>


template <typename Func, typename... Args>
class CallH {
public:
    // Helper function to call any function
    static void call(Func&& func, Args&&... args){
    std::forward<Func>(func)(std::forward<Args>(args)...);
    };
};

#endif // CALLABLE_HELPER_H

#pragma once
#include <string>
#include <iostream>

#define _MAKESTR(x) #x
#define MAKESTR(x) _MAKESTR(x)

inline void msg_assert(const char* condition, const char* message, const char* fileline) {
    std::cerr << "[" << fileline << "] "
        << "Assertion `" << condition << "` failed.\n"
        << message << std::endl;
    std::abort();
}

inline void msg_assert(const char* condition,
    const std::string& message,
    const char* fileline) {
    msg_assert(condition, message.c_str(), fileline);
}

#ifdef NDEBUG
#define better_assert(condition, message) static_cast<void>(0)
#else
#define RgAssert(condition, message) \
    static_cast<bool>(condition)          \
        ? static_cast<void>(0)            \
        : msg_assert(#condition, message, __FILE__ ":" MAKESTR(__LINE__))
#endif
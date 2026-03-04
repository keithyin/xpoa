cmake_minimum_required(VERSION 3.10)
project(xpoa VERSION 1.0 LANGUAGES CXX)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. 搜集所有源文件
# 包含根目录下的 .cpp 和 edlib 目录下的 .cpp
file(GLOB_RECURSE ALL_SOURCES 
    "*.cpp"
    "edlib/*.cpp"
)

# 2. 区分库文件和可执行文件源文件
# 从总列表中移除 main.cpp，剩余的用于构建静态库
set(LIB_SOURCES ${ALL_SOURCES})
list(FILTER LIB_SOURCES EXCLUDE REGEX "main\\.cpp$")

# 3. 创建静态链接库
# 库名设为 xpoa_lib
add_library(xpoa STATIC ${LIB_SOURCES})

# 指定头文件包含路径，确保能找到根目录和 edlib 目录下的 .h
target_include_directories(xpoa PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/edlib
)

# 4. 创建可执行二进制文件
add_executable(xpoa_bin main.cpp)

# 5. 将静态库链接到可执行文件
target_link_libraries(xpoa_bin PRIVATE xpoa)
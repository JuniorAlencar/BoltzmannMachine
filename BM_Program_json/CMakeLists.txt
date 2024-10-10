cmake_minimum_required(VERSION 3.10)

# Nome do projeto e versão do C++
project(MeuProjeto VERSION 1.0)

# Usar Clang como compilador
set(CMAKE_CXX_COMPILER "clang++")

# Definir flags de compilação
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -D_GLIBCXX_USE_CXX11_ABI=1 -Wno-deprecated -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_BUILD_TYPE Debug)

# Adicionar flag para desabilitar a inclusão de cabeçalhos desatualizados
add_definitions(-D_LIBCPP_ENABLE_CXX17_REMOVED_AUTO_PTR)

# Adicionar todos os arquivos .cpp, .hpp e .h do diretório src/
file(GLOB_RECURSE SOURCES src/*.cpp src/*.hpp src/*.h)

# Incluir bibliotecas Boost (filesystem e system)
find_package(Boost REQUIRED COMPONENTS filesystem system)

# Verificar se o Boost foi encontrado
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(meu_executavel ${SOURCES})
    target_link_libraries(meu_executavel ${Boost_LIBRARIES})
else()
    message(FATAL_ERROR "Boost not found")
endif()


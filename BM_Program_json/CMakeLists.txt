cmake_minimum_required(VERSION 3.10)

# Nome do projeto e versão do C++
project(MeuProjeto VERSION 1.0)

set(CMAKE_CXX_COMPILER "g++")

# Set the C++ 14 standard
set(CMAKE_CXX_STANDARD 14)
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Add -O3 optimization flag
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Incluir diretórios
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Adicionar arquivos fonte (excluindo main.cpp)
file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.hpp" "src/*.h")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

# Adicionar o executável para a aplicação principal
add_executable(meu_executavel src/main.cpp ${SOURCES})

# Incluir bibliotecas Boost (filesystem e system)
find_package(Boost REQUIRED COMPONENTS filesystem system)

# Verificar se o Boost foi encontrado e imprimir a versão
if(Boost_FOUND)
    message(STATUS "Using Boost version: ${Boost_VERSION}")
    target_link_libraries(meu_executavel PRIVATE ${Boost_LIBRARIES})
else()
    message(FATAL_ERROR "Boost not found")
endif()

# Set the output directory for the main executable
set_target_properties(meu_executavel PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# Desativar o uso de auto_ptr no Boost
add_definitions(-DBOOST_NO_AUTO_PTR)

# Forçar a ABI correta do GCC para evitar uso de auto_ptr
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)

# Encontrar todos os arquivos de teste
file(GLOB TEST_FILES "tests/*.cpp")

# Criar um executável separado para cada arquivo de teste
foreach(TEST_FILE ${TEST_FILES})
    # Obter o nome do arquivo sem a extensão
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    
    # Adicionar o executável para cada teste
    add_executable(${TEST_NAME} ${TEST_FILE} ${SOURCES})

    # Linkar as bibliotecas necessárias
    target_link_libraries(${TEST_NAME} PRIVATE ${Boost_LIBRARIES} pthread)

    # Definir o diretório de saída para os executáveis de teste
    set_target_properties(${TEST_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    )

    # Registrar o teste
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()

# Habilitar testes
enable_testing()

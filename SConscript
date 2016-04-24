# SConscript

Import('env Automatic_Program')
env=env.Copy(warnings_are_errors=1)
Automatic_Program(env,"embedded_elasticity",
                  [
                      "main.cpp",
                      "ELASTICITY_DRIVER.cpp",
                      "ELASTICITY_EXAMPLE.cpp",
#                      "ELASTICITY_EXAMPLE_2D.cpp",
                      "ELASTICITY_EXAMPLE_3D.cpp",
                      "EMBEDDED_SURFACE.cpp",
                      "EMBEDDED_SURFACE_3D.cpp",
#                      "EMBEDDED_SURFACE_2D.cpp",
#                      "EMBEDDING_TOOLS_2D.cpp",
                      "EMBEDDING_TOOLS_3D.cpp",


                  ])


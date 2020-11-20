# Tareas para los cambios en el paquete

## Carga de datos _single-cell_ en el objeto `DigitalDLSorter`

1. Hacer una función para convertir CSV demasiado grandes en ficheros HDF5. Mirar esta página de stackoverflow donde lo hacen en Python, asumo que se puede hacer de la misma manera en R:
    * <https://stackoverflow.com/questions/27203161/convert-large-csv-to-hdf5>
    * <https://stackoverrun.com/es/q/7467814>

    La función podría trabajar sí o sí independientemente del tamaño de los datos o dar un argumento al usuario para poderlo usar? No sé. Hay que tener en cuenta que lo lógico es que los datos sean cargados en el objeto a través de la función `loadSCProfiles` sí o sí.

2. Parámetro `file.backend` donde se localizará el fichero HDF5 con los datos.

## Carga de datos _single-cell_ en el objeto `DigitalDLSorter` (versión definitiva)

* La función va a permitir cargar matrices de expresión en formatos TSV, CSV, matrices sparse y ficheros HDF5.
* El hecho de utilizar o no ficheros HDF5 como backend depende el argumento `back.end`. Aunque se introduzcan ficheros HDF5 como input, solo se usarán como backend en el caso de que se dé en ese argumento un path válido. 
* En el caso de que no se utilicen ficheros HDF5 pero se quieran utilizar como _back-end_, bastará con proporcionar un path correcto en el argumento `file.backend`.
* Finalmente, en el caso de que se quiera trabajar por bloques (solo compatible con el uso de ficheros HDF5 como back-end), habrá que especificar el argumento `block.processing = TRUE`.

## Función `estimateZinbParams`

* Argumento para decidir el número de células que van a formar el subset (tiene que ser menor que el número de células totales)
* Argumento para decidir si el número de células de cada tipo va a ser igual o si se van a guardar las proporciones originales. 
    * En cualquier caso, si hay algún tipo celular que no llegue al número deseado, avisar de alguna manera, pero no devolver un error.

## Función `simSCProfiles` (cambiar nombre)

* Falta permitir la posibilidad de utilizar ficheros HDF5 como backend (fácil).
* Falta la implementación del procesamiento por bloques (complicado). Respecto a esto último, la idea es hacerlo con un bucle for y utilizar el  paquete rhdf5 que, con suerte, permitirá la adición de nuevos elementos en un fichero HDF5 preexistente. Se podría intentar implementar en C++, pero habría que mirar un paquete estadístico que implemente distribuciones binomiales. 
    * <https://support.bioconductor.org/p/69373/>: explican cómo añadir filas a un dataset en un fichero HDF5.
    * En el caso de implementarlo en C++ (no merece la pena), hay gente con el mismo problema: <https://stackoverflow.com/questions/15379399/writing-appending-arrays-of-float-to-the-only-dataset-in-hdf5-file-in-c>, <https://forum.hdfgroup.org/t/c-interface-append-data-to-an-existing-file/4474>.

## Función `generateBulkCellMatrix` (cambiar nombre)

Modificaciones: 

* Introducir la posibilidad de establecer la proporción de muestras bulk de train y de test (ahora el argumento `train.freq` solo afecta a las células).
* Eliminar lo de los tipos celulares exclusivos.
* Fusionar tablas de metadatos.
* A la hora de dividir las células en train y test, se hace con el conjunto completo de células. Dado que el número inicial de células puede ser reducido, quizás sea conveniente introducir un argumento rollo `balanced = TRUE` para dividir las células en train y test por clases. Voy a meterlo, pero, en cualquier caso, no puede fallar como lo hace actualmente. 

## Desarrollo del entrenamiento de la red neuronal 'on the fly'

Cosas a tener en cuenta:

* La idea es que los generadores llamen a las funciones necesarias para simular las muestras _bulk_.
* Hacer shuffling sobre el nombre de las muestras en lugar de usar la matriz.
* Permitir modificar la arquitectura de la red: permitir modificar el número de capas y el número de neuronas mediante los argumentos y permitir añadir un modelo compilado con la arquitectura deseada (siempre y cuando sea una red neuronal completamente conectada).

## Dudas durante el desarrollo

* Qué hacer con el argumento `setType` de la función `estimateZinbParams`: es para establecer un tipo celular a evaluar. **Implementado**.

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


## Dudas durante el desarrollo

* Qué hacer con el argumento `setType`de la función `estimateZinbParams`: es para establecer un tipo celular a evaluar.

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

## Función `generateBulkCellMatrix` (cambiar nombre) --> simulBulkSamples

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
* Probar diferentes versiones para llevar a cabo la operación de rowSums. Algunas ideas usando `dplyr` o `data.table`: <https://stackoverflow.com/questions/15905257/getting-rowsums-in-a-data-table-in-r>, <http://brooksandrew.github.io/simpleblog/articles/advanced-data-table/>.


**Novedades**

Está implementado. Las siguientes cosas que hay que implementar son:

* Probar un apply más rápido para acelerar el proceso (o incluso implementarlo en C++, igual es algo que debemos hacer).
* Hacer una función que permita mantener las muestras bulk en ficheros.
* Probar otras arquitecturas y probar meter un modelo keras ya hecho.
* Hay que modificar los slots de la clase para añadir las frecuencias target: dado que lo haces on the fly, necesitas hacerlo en algún momento.
* Quitar el cálculo de las métricas por tiempo e intentar hacerlo a mano o mirar alguna alternativa que dé Keras teniendo las matrices. Si no, pues nada.

**Importante**

Las métricas de evaluación se pueden calcular a mano para reducir los tiempos de ejecución en vez de hacerlo dos veces. Mirar a ver si renta, porque veo varios inconvenientes respecto a manejar cómo coger la función correcta cuando el usuario mete las métricas de evaluación en el argumento `metrics`.

## Función `simulBulkProfiles`

* Utilizar la misma estructura que a la hora de simular los perfiles single cell: dar la opción de hacerlo en memoria y mantenerlo en memoria, hacerlo en memoria y pasarlo a un fichero HDF5 y hacerlo por bloques en un fichero HDF5. 
* Basta con utilizar la función `setBulk` según convenga: o por chunks o sin ellos. Dado que voy a hacer el shuffling antes de generar las muestras independientemente de cómo se vayan a construir, el slot `bulk.simul` debería contener el orden de éstos (al menos en el entrenamiento). Cuando se hace _on the fly_ realmente da lo mismo, porque esas muestras no las va a volver a ver, solo se guarda el de test y no hace falta porque no se hace shuffling en el test.
* Habrá que construir versiones nuevas de los generadores para entrenar a la red neuronal o modificar los ya hechos. Creo que puede rentar modificarlos, únicamente bastaría con una sentencia `if` con un parámetro indicando si construir las muestras correspondientes a los perfiles de la variable `sel.data` o cogerlas del slot correspondiente. El tema de pillar los perfiles single-cell de los slot correspondientes ya está implementado correctamente.

**Importante**

Dado que solo voy a implementar la posibilidad de generar las matrices _bulk_ sin la normalización-transposición antes de entrenar a la red neuronal, sería conveniente que dichos pasos se hagan de forma independiente en los generadores de matrices. De esta forma, al ser tareas independientes y hacerse de forma modular, queda mejor.

**Problema**

* Como me temía, el hecho de tener las células simuladas en disco ralentiza enormemente el proceso hasta el punto de tardar para 500 muestras unos 17 minutos. Voy a utilizar algunos de los parámetros del paquete `rhdf5` para intentar acelerar el proceso, pero no sé hasta qué punto funcionará.
* Optimizando el argumento `chunk`, para 500 muestras pasamos de 17 minutos a 11. El siguiente paso va a ser ordenar la forma de acceder.
* Ordenando tampoco mejora. 1000 muestras son 20 minutos. El siguiente paso es cambiar el tamaño de chunk y probar accediendo mediante las funciones del paquete `rhdf5`.

**Solución**

Finalmente 1000 muestras tardan alrededor de 5 minutos ordenando las células y utilizando un chunk = c(27mil, 1). No sé si habrá una combinación más eficiente.

## Función `trainDigitalDLSorterModel`

Dos versiones por debajo: una para el entrenamiento _on the fly_ y otra si las muestras bulk están ya generadas.    

## Dudas durante el desarrollo

* Qué hacer con el argumento `setType` de la función `estimateZinbParams`: es para establecer un tipo celular a evaluar. **Implementado**.
* Mirar si el hecho de cargar entero un paquete dentro del paquete produce lentitud. En tal caso, `keras` probablemente sea un paquete relativamente grande, por lo que será mejor llamar a las funciones explícitamente.
* Hay muestras Bulk llenas de ceros. Esto es un problema a la hora de generarlos. Generar todo de nuevo y ver si el problema persiste.


## Links de interés

* Explican lo del argumento `chunk`: <https://www.bioconductor.org/packages/devel/bioc/vignettes/rhdf5/inst/doc/practical_tips.html>
* Link para implementar lo de meter filas en la matriz: <https://stackoverflow.com/questions/37625471/insert-row-at-a-specific-location-in-matrix-using-r>

## Tareas que quedan

* Probar todas las nuevas funcionalidades: implementar el preprocesamiento por chunks (mirar los libros de los ficheros HDF5).
* Probar diferentes chunk.dim para comprobar cuál es la combinación más rápida en las muestras bulk.
* Mirar lo de simular los perfiles de expresión con regresión linear.
* Mirar el cálculo de los saliency.

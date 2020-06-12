# Installation
Read manual

# Compiling
To compile the program simply do:

    make

This will compile the program and move the executable to the 'bin' folder.

To compile and run simply do:

    ./script
    
# Selección de espacio de color, región y filtros
Para seleccionar el punto final de la ruta se da click derecho en la imagen.

Para seleccionar el espacio de color, la imagen debe estar pausada. Para pausar la imagen se da en la tecla de espacio. Una vez pausada se puede con el mouse dibujar un rectángulo con el espacio de color.

Una vez seleccionado el espacio de color con la tecla f se filtra. Por defecto se filtra usando BGR pero con la tecla y se convierte a YIQ y con la tecla v a HSV.

# Entrenamiento
Para el entrenamiento se tiene que descomentar la línea de código con el #define TRAIN

Una vez binarizada la imagen al presionar r se buscan las regiones y se calculan los momentos de Hu.

Al presionar n, se guardan las phiś de ese instante como muestra y se pasa a la siguiente muestra.

Esta programado para guardar 10 muestras por objeto y son 4 objetos. Al presionar n por décima vez, se cierran las ventanas y se escribe el promedio y la desviación de las phi's del objeto en el archivo "parameters.txt". Se repite el proceso para los 4 objetos.

# Corriendo el clasificador
Después de la selección del espacio de color, región y filtros, con la tecla r se calculan las phi's y se comparan con parameters.txt. En la ventana de la mira se prende el cuadrante correspondiente y la flecha sigue la dirección del objeto largo.

Con la tecla k se calcula el camino y se muestra en pantalla.

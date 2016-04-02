# Projet de bio-informatique

## Structure du projet

```
README.md     // ce fichier
build.xml     // fichier pour ant
src           // dossier contenant les fichiers sources
|-->  main    // sources exceptés les tests.
|-->  test    // sources des tests.

bin           // dossier contenant les .class. Généré par ant cmpile.

lib           // dossier contenant les librairies externes comme junit.
```

Chaque package source commence par *be.ac.umons.bioinfo* et les tests également.

Les fichiers sources se trouvent dans le dossier *src/main* tandis que les fichiers tests dans *src/test*.

Le nom de package du fichier *src/test/be/ac/umons/bioinfo/PACKAGE/ATest.java* est le même que *src/main/be/ac/umons/bioinfo/PACKAGE/A.java* qui est *be.ac.umons.bioinfo.PACKAGE*.

## Utilisation du build.xml avec Ant

* **ant**: lance *compile*
* **ant compile**: compile les fichiers sources du projet, excepté les fichiers
  tests.
* **ant compile_test**: compile les fichiers tests qui se trouvent dans le
  dossier *src/test* après avoir lancé la cible *compile*.
* **ant run_test**: lance les tests.
* **ant run**: lance le main après avoir lancé run_test.
* **ant zip**: crée une archive zip avec les fichiers sources, le build.xml et le README. Sont exclus les dossiers et fichiers .git, bin, .project, .classpath, .settings/ et todo.
* **ant clean**: supprime le dossier bin généré lors de la compilation.
* **ant fclean**: lance la dépendance *clean* et supprime l'archive *zip* si créée.

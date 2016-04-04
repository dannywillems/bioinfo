# Projet de bio-informatique

[![Build Status](https://travis-ci.org/dannywillems/bioinfo.svg?branch=master)](https://travis-ci.org/dannywillems/bioinfo)

Testé avec:
* Open JDK 6
* Open JDK 7
* Oracle JDK 7
* Oracle JDK 8

## Structure du projet

```
README.md             // ce fichier
build.xml             // fichier pour ant
src                   // dossier contenant les fichiers sources
|-->  main            // sources exceptés les tests.
|-->  test            // sources des tests.

bin                   // dossier contenant les .class. Généré par ant cmpile.

lib                   // dossier contenant les librairies externes comme junit.

doc                   // dossier contenant les fichiers de documentations. *ant doc*

res                   // ressources du projet
|--> collections      // collections de fichiers fasta
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
* **ant zip**: crée une archive zip avec les fichiers sources, le build.xml et le README. Sont exclus les dossiers et fichiers .git, bin, .project, .classpath, .settings/, .idea, .travis.yml, bioinfo.iml.
* **ant doc**: construit la documentation grâce à javadoc.
* **ant clean**: supprime le dossier bin généré lors de la compilation.
* **ant fclean**: lance la dépendance *clean* et supprime l'archive *zip* si créée.

## To-do

- [ ] Tests unitaires sur lecture/écriture de fichier fasta.
- [x] Lectures et écritures de fichier fasta

- [x] Implémentation de l'alignement semi-global
- [x] Tests unitaires de l'alignement semi-global

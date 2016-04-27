# Projet de bio-informatique

[![Build Status](https://travis-ci.com/dannywillems/bioinfo.svg?token=VX8gT1NE5x87pjz7pBvN&branch=master)](https://travis-ci.com/dannywillems/bioinfo)

Testé avec:
* Oracle JDK 8

**Nécessite JDK >= 8.**

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

#### Fichiers fasta
- [ ] Tests unitaires sur lecture/écriture de fichier fasta.
- [x] Lectures et écritures de fichier fasta.

#### Alignement
- [x] Implémentation de l'alignement semi-global.
- [x] Tests unitaires de l'alignement semi-global.

#### Algorithme greedy
- [x] Implémenter l'algorithme.
- [ ] Multi threader l'algorithme.
- [ ] Fonctionne quand on a une séquence qui est dans l'autre ?
- [ ] Unit test pour l'algorithme greedy.

#### MISC
- [ ] Suppression des tests *inutiles*.
- [ ] Suppression des commentaires *inutiles*.

## Plus d'informations

* Introduction to computational molecular biology - Carlos Setubal, Joao Meidanis. Chapitres 3 et 4.

#### Alignement de séquences

* Wikipedia: https://en.wikipedia.org/wiki/Sequence_alignment

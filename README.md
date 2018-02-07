# Xylo

Arficial life simulation of plant-like cells.

**TL;DR for nerds:** Xylo evolves cells with GRN-like dynamics and
epigenetic marks for tissue specialisation. Cells use fuzzy string
matching (like [StringMol](http://stringmol.york.ac.uk/)) to react to
stimuli and internal products.

The simulation takes place in a simple world, currently 2D. Cells
occupy a position in space and, as you would expect, are subject to
gravity and collide with other cells and the ground. They may bond to
adjacent cells, thus forming multi-cellular structures, and may break
such bonds. They may exert a force (accelerate) in a direction. They
may reproduce, either asexually or sexually. All these operations
consume energy.

Being plant-like cells, they generate energy from sunlight, meaning
that they continuously generate energy unless shaded from
above. Actually let's call it sugar, since cells can accumulate it up
to a certain amount. It is also possible for a cell to open a sugar
channel to an adjacent cell in order to pass energy, much as the
leaves of a tree pass sugars down to its roots.


## Computation

Despite the term _cell_, these units are in some ways more like the
molecules of artificial chemistries.  A cell is characterised
primarily by its "DNA", a sequence of bases. Unlike most artificial
life models, this DNA has epigenetic properties: parts can be
silenced. The non-silenced parts are called _open_ DNA. When a cell
comes in contact with an external stimulus, the stimulus is always
presented as a DNA-like sequence which is matched against the cell's
open DNA. In addition, cells allow internal dynamics by accumulating
molecular products -- more DNA-like sequences -- which again match
against the cell's open DNA. Internal dynamics just means that there
is some memory; behaviour can be more complex than an immediate
reaction.

A match to a section of DNA can be inexact: a score is calculated
allowing for mismatches, insertions and deletions (the Smith-Waterman
algorithm). Since there may be multiple possible matching sites across
a DNA sequence, these similarity scores are used to weight the
probability of selecting each site. Once a site is bound, the DNA is
read from that point forward and interpreted as a sequence of
operations. These are introduced below.

### Reproduction

a/sexually

tissue specialisation
silencing


### Physical operations

* apply force
* form/break bond
* about face
* rotate left/right

### Sugar

* start/stop a sugar channel

### Epigenetics

* silence/unsilence DNA

### Regulatory networks

* create molecular product

- products decay


### Control flow

* stop this reaction
* jump to another DNA location
* conditional check energy level




## Death




## Stimuli

* ground
* sunlight?
* other cells




## Usage

FIXME

## License

Copyright © 2018 Felix Andrews

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.

# issues

viz products

viz sugar channels

enforce spatial boundaries

push into cell could destroy other cell? - attack & consume?


## impl


* check for cell in current orientation as goto-like


hand-craft genomes to test each operation.


## design

better for each stimulus to only bind to the single best match -- predictable masking.

adjust match indexes by 1 (possible match to start of sequence)

internal dynamics (GRN) should have negatives -- repressors -- as well as positive TFs


what is the main loop for a cell?
  - reacting to new external stimulus only?
    - though including results of our own actions: reproduction / forces / sugar channels
  - what about internal dynamics?
    - with GRN, input/output seemed awkward
    - with stringmol matching, i/o is natural but no internal dynamics.
    - delay-goto - carry over a reaction to next time step
      - but vs new stimuli? and what about infinite loops.
    - product - accumulates *inside* cell and *competes* to react.
      - also decays.
      - produce up to a max level of each, and have a max number of different TFs.
      - i.e. GRN-like.
      - then epigenetic silencing could be limited to (and mandatory for) reproduction
        - tissue specialisation

*********
alignment should be stochastic based on match score
and considering all candidate alignment locations.
*********

matching should include a monotonic function of alphabetical distance per base
  - as a continuous selection signal
  - alternatively: reduce alphabet to ~4 bases and interpret codons as operations.
    - correspondingly increase match score threshold.

_physically oriented operations_
1. 0 to-stimulus
1. 0 to-sun
1. 0 about-face
1. 0 rot-left
1. 0 rot-right
1. 4 clone
1. 4 sex - when stimulus is another cell
1. 2 push
1. 1 bond-form   }
1. 1 bond-break  } when there is an adjacent cell in current orientation
1. 1 sugar-start }
1. 1 sugar-stop  }

_non oriented operations_
1. 1 product - read template (up to next *terminator* or *stop*) - min length is 2 codons.
1. 0 silence - generates a product which silences where it binds - from start to end of match.
1. 0 unsilence - generates a product which unsilences where it binds _against full DNA_
1. 0 goto - read template, find matching site as like a product (with translation), transfer control.
1. 0 energy-test - threshold is proportional to current read location; same as goto _if_ energy meets threshold.
1. 0 terminator (a set of codons?)
1. 0 stop reaction

= 18 ops + terminator(s).

codons of 3 in alphabet of 4 = 64.
  - leaves 46
    - if 12 are terminators, expected length of templates is ~5 codons.
  - distributions throughout encoding space?
    - need non-terminator templates to match to terminators - in both complement and (product) translation.
      - maybe not necessary / can always match approximately

is binding/unbinding once per stimulus?
  - if created touching rock, only react once?
  - or do we re-evaluate if a conditional changes / or epigenetic state changes?
    - fine if the reactions are idempotent
      - but they are not: reproduction / push / silence

all reactions should by default have a persistent effect

so, cell consists of genome + epigenetic marks + TFs
  - stimuli and TFs compete to react each step
    - with probability according to goodness of match.
  - TFs decay at some rate

how does unbinding react anyway?
  - reverse sequence for alignment?

once a sugar channel is opened, does it
  - pass sugar every step (until threshold?)
  - react every step to query

energy-test conditional
  - threshold is calculated from position in genome from 0 to max energy level

limit instructions per step
  - at most one copy per step.
  - children can not copy in first step. (anyway: energy cost)

energy cost of reactions - esp clone


matching to self clone:
  - match to complement is symmetric so two equally good candidate binds on each side
  - use first in sequence to look up reactions
  - need specialisation to be interesting
    - is silencing enough?
    - execution forks at copy time.
      - run through parent first. when that ends, child begins?
      - or, child could be anti-sense?
    - want/expect silencing in parent and/or child.
      - make this the default behaviour - always silence by reading next templates
        - read begin token up to next opcode (any). then end token up to following opcode.
        - child uses following two templates, similarly.
      - unsilencing?
        - a normal action
    - additional specialisation mechanism:
      - initialise new cell with special birth products


also: GRN-like products are transcribed from a template in DNA
  - product will just match its original template
  - so need translation
  - if in complement, the template will match binding site in a clone cell
  - so need translation different to cell-cell matching (not complement)

silencing
  - continue to best match template alignment
  - wrap around?
    - or only consider downstream genome
    - or use current direction (fwd/back)
  - if no match (0), no-op

unsilencing
  - match to silenced bases of genome?

sexual reproduction - align genomes at crossover points (use trace)

## future

water

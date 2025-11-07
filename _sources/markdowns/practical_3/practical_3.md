## Setting conservation priorities

### 1. Introduction and resources

This practical is aimed to introduce you to EDGE scores that you’ll need
for your conservation strategy coursework. Put briefly, these scores
balance the distinctiveness of species against their risk of extinction
to determine conservation priorities. You can find out more information
about EDGE scores from the ZSL website:

<https://www.zsl.org/conservation/our-priorities/wildlife-back-from-the-brink/animals-on-the-edge>

### 2. Preparing data

To calculate EDGE metrics, we need data on the species we’re interested
in, and their phylogenetic relationship. For the coursework we’re
interested in EDGE scores for a specific country, however it’s also
common to look at specific clades.

For this practical we’re going to use data from Iceland. To do this we
will need a list of bird species in Iceland and the Avonet data used in
practical 2.

Because the 2 data sets may not line up we will also load a crosswalk
data set to easily merge our focus data sets.

``` r
library(phytools)
library(phylobase)
library(dplyr)
library(stringr)
library(ape)

#import all bird species
avonet_data <- read.csv("data/avonet_data.csv")

#import country species conservation data
conservation_data <- read.csv("data/Iceland-Species.csv")

#import crosswalk 
crosswalk <- read.csv('data/BL_Jetz crosswalk.csv')
```

Use the crosswalk to prepare the Avonet data filter.

We should standardise the species names across all of our data sets.
Some vary between having a space between them and an underscore, this
will set them all to an underscore

``` r
conservation_data$ScientificName <- str_replace_all(conservation_data$ScientificName, " ", "_")
crosswalk$BirdLife.name = str_replace_all(crosswalk$BirdLife.name, ' ', '_')
avonet_data$birdlife_name = str_replace_all(avonet_data$birdlife_name, ' ', '_')
```

We can now filter the Avonet data for each species found in your
selected country, this will give us trait data and redlist status of
each species which will be useful when calculating extinction
probabilities.

``` r
#filter for only country species
Species <- avonet_data %>%
  filter(birdlife_name %in% conservation_data$ScientificName)
```

We now need to clean the data to avoid any possible errors with tree
formation later on.

First we need to remove any species that are extinct or extinct in the
wild before we calculate extinction probability. These species would
show with the highest EDGE score as they would have an extinction
probability of 1.

Second we should remove all species that are data deficient or NA, this
will cause our tree and data to misalign which will cause our ancestral
calculations to be incorrect.

In addition to this depending on your selected country you may need to
remove all Procellariiformes, these are usually included in most
countries due to their large ranges but do not land or nest there.

It’s important to remove all unwanted data before creating our tree as
this will cause mismatches later but won’t create an error message. For
this reason it’s also important to check that your results make sense,
if your tree is misaligned your results could change each time you run
the code.

``` r
#clean the data 

# - remove extinct species, if needed
Species <- Species %>% 
  filter(!(redlist_cat %in% c('EX', 'EW', 'DD')))

# - remove NAs, if needed
Species <- na.omit(Species)


#remove sea birds if they dont nest/ stop in selected country
#Species <- Species %>% 
 # filter(jetz_order != 'Procellariiformes')
```

### 3. ED Scores

Instead of evolutionary distinctiveness, which would be conventional in
EDGE2 calculations, we might instead be interested in what functions
each species provides the ecosystem. Species with low functional
distinctiveness may be ‘functionally redundant’ in the ecosystem,
whereas those with high functional distinctiveness may provide key
ecosystem services that aren’t easily replaceable.

To calculate the function of species, we’ll use their morphological
traits. Recent research (Pigot *et al.* 2020) has shown that a few
simple traits can predict the niche a species occupies, which in turn
tells us about function. For instance, frugivores provide vital seed
dispersal, and insectivores influence invertebrate population dynamics.

You can read the full article here:
<https://www.nature.com/articles/s41559-019-1070-4>

We’ll use just a few key traits, separated into beak traits (primarily
trophic), body traits (primarily locomotion), hand wing index (flight
ability), and body mass.

First we should check correlations between our traits:

``` r
# First log body mass.
Species$logMass <- log(Species$mass)

# Split the traits into beak shape traits, body shape traits, and body mass.
beak_traits <- Species %>% dplyr::select(beak_length_culmen, beak_width, 
                                            beak_depth, logMass, hand_wing_index)
body_traits <- Species %>% dplyr::select(tarsus_length, wing_length, 
                                            tail_length, logMass, hand_wing_index)

# Look at the correlations between traits, including body mass.
pairs(beak_traits)
```

```{image} practical_3_files/figure-gfm/unnamed-chunk-6-1.png
:align: center
:width: 600px
```

``` r
pairs(body_traits)
```

```{image} practical_3_files/figure-gfm/unnamed-chunk-6-2.png
:align: center
:width: 600px
```

So we can see some pretty strong correlations between beak shape, body
shape and body mass. This is because bigger species tend to have
proportionally bigger traits. If we were to calculate functional
distinctiveness using these traits in their current form, we would
probably only determine which species are biggest and smallest. Another
problem is that when all the traits are highly correlated, it’s harder
to tell which ones are important for functional distinctiveness. Instead
our model may simply pick them at random, leading to potentially random
results.

Instead we can use two variables from Pigot *et al.* (2020), “beak
shape”, and “body shape”. These variables were created using a Principal
Components Analysis (PCA) to condense our traits down into fewer
uncorrelated variables, removing body mass as the dominant effect.

```{tip}
We won't go into much detail on how a PCA works, as we want to focus on understanding the role of functional distinctivness in conservation. If you'd like to know in more detail, StatQuest has a great 5 minute [video](https://www.youtube.com/watch?v=HMOI_lkzW08&ab_channel=StatQuestwithJoshStarmer).
```

After removing body mass as the dominant effect, beak shape correlates
positively with beak length, and negatively with beak width and depth.
So a higher value beak shape means a narrower longer beak. For body
shape, there is a positive correlation with tarsus length, and a
negative correlation with tail length. So the most extreme would be tall
birds like Ostrichs and Flamingos.

We can check to see if there’s still any strong positive or negative
correlations:

``` r
#create trait dendrogram

functional_traits <- Species %>% 
  dplyr::select(beak_shape,body_shape,logMass,hand_wing_index)

# Look at the correlations between traits.
pairs(functional_traits)
```

```{image} practical_3_files/figure-gfm/unnamed-chunk-8-1.png
:align: center
:width: 600px
```

No clear positive or negative patterns, so we’re okay to go forward and
calculate functional diversity.

First we need to create a distance matrix of our traits. Species with
similar traits will have smaller ‘distances’.

``` r
#create distance matrix based on functional traits
functional_traits <- functional_traits %>% 
  scale()
rownames(functional_traits) <- Species$jetz_name
dist_matrix <- dist(functional_traits)
```

The next step is to create a tree using the neighbour-joining method
(Saitou & Nei, 1987) (Google for more information!). This will create a
tree where branch lengths show how similar species are in trait space
rather than evolutionary distance. This function may take a while with
more species so don’t be alarmed if the country you’ve chosen takes much
longer.

``` r
#build a neighbour-joining tree based on trait distances
traittree <- nj(dist_matrix)

#root the tree as its midpoint
tree_nj <- midpoint_root(traittree)

#plot the tree

plot(traittree, cex=0.3)
```

```{image} practical_3_files/figure-gfm/unnamed-chunk-10-1.png
:align: center
:width: 600px
```

```{tip}
Have a play around with the "plot.phylo()" function in the {ape} package. Try and see what alternative ways you can present this tree.
```

### 4. Extinction Probability

Now that we have a tree for our ED we can move onto setting our
extinctinion probabilities.

More complex calculations, such as in (Gumbs et al,2023), assign a range
of extinction probabilities to each redlist category to cover the
uncertainty of extinction.

To do this they create 1000 pext values for each GE and randomly select
one for each species, this means that they have to iterate through the
code until they get the most likely EDGE scores. While this is more
accurate it would take significantly longer than our code.

You can read the full paper here if you’re interested:
<https://journals.plos.org/Plosbiology/article?id=10.1371/journal.pbio.3001991>

``` r
#generate extinction probabilities for each redlist cat - if error check Species$redlist_cat 
pext.vals <- data.frame(rl.cat = rev(c("CR","EN","VU","NT","LC")) , 
                        pext = rev(c(0.97, 0.97/2, 0.97/4,0.97/8,0.97/16)))
Species$pext <- pext.vals$pext[match(Species$redlist_cat, pext.vals$rl.cat)]

Refinedspecies <- data.frame(species = Species$jetz_name, pext = Species$pext)

head(Refinedspecies)
```

    ##                species     pext
    ## 1 Haliaeetus_albicilla 0.060625
    ## 2           Anas_acuta 0.060625
    ## 3        Anas_clypeata 0.060625
    ## 4          Anas_crecca 0.060625
    ## 5        Anas_penelope 0.060625
    ## 6   Anas_platyrhynchos 0.060625

### 5. EDGE2 Scores

To make future necessary calculations more robust we should write a
function to order our tree when needed.

``` r
# order tree components
reorder_tree <- function(tree, ordering){
  tree@edge.length <- tree@edge.length[ordering]
  tree@edge <- tree@edge[ordering,]
  return(tree)
}
```

We need to change the format of our tree from “phylo” to “phylo4”, this
allows us to use better functions to manipulate our tree.

We also need to define how many species we have, how many nodes we have,
and how many data points in total we have. This will allow us to fully
calculate the tree length for each species including it’s ancestors
length.

``` r
tree <- as(tree_nj,'phylo4')
root <- rootNode(tree)
N_species <- length(tree_nj$tip.label)
N_nodes <- tree_nj$Nnode
N_tot <- N_species + N_nodes
```

Once we’ve done this we can create a dataframe to store our results

``` r
#create a data frame to store all calculated data
tree_data <- data.frame(species = as.character(tree_nj$tip.label),
                        TBL = tree@edge.length[1:N_species],
                        pext = Refinedspecies$pext, 
                        ED = NA, 
                        EDGE = NA)
```

We can use the function we just wrote to align our tree and nodes, this
will mean that our ancestral tree length will always be calculated
before our species length. If we skip this step we will end up with much
shorter tree edge lengths as the species will have an ancestral length
of 0.

``` r
nodes <- c(root,descendants(tree, root, 'all'))

ord <- order(nodes)
tree <- reorder_tree(tree, ord)
nodes <- nodes[ord]
#calculate the edge length of each species in our tree
tree_data$TBL <- tree@edge.length[1:N_species]

node_data <- data.frame(Node = 1:N_tot, Pext = as.numeric(rep(1, N_tot)), Edge_Sum = NA)
node_data[1:N_species, 2] <- as.numeric(Species$pext)
```

Now that everything is ordered we can begin to calculate the EDGE score
for our species.

``` r
# assign the product of its descendant tips to each node
for (i in c(1:length(tree@label), N_tot:(root+1))){      # for each node, beginning with tips
  ans <-tree@edge[i,1]                                       # find ancestor of node
  node_data[ans, 2] <- node_data[ans, 2]*node_data[i,2]     # muliply ancestor value by node "pext"
}

# multiply each edge.length by each pext calculated above
for (i in 1:length(nodes)){
  tree@edge.length[i] <- tree@edge.length[i]*node_data[i,2]
}


if (is.na(tree@edge.length[root])){
  tree@edge.length[root] <- 0
}
```

When we created the tree we rooted at the midpoint, this is because the
tree needs to be rooted for the following calculations and we do not
have an outgroup to root the tree. We also need to set the root length
to 0 as it is at the base of the tree and then will allow us to measure
the length of all of our species.

In our tree we have a list of each ancestor and their descendants, using
this and our root we can summate the ancestral edge length of our
species. By nesting this within a while loop we can control how long the
loop continues, Once we reach the bottom of the list the loop
terminates.

``` r
# for each internal node, summate ancesteral edgelengths
node_data$Edge_Sum[root] <- 0

change <- TRUE
while(change) {
  change = FALSE
  for( i in seq_len(nrow(tree@edge))) { 
    parent = tree@edge[i,1]
    child = tree@edge[i,2]
    if(!is.na(node_data$Edge_Sum[parent]) && is.na(node_data$Edge_Sum[child])) {
      node_data$Edge_Sum[child] = node_data$Edge_Sum[parent] + tree@edge.length[i]
      change = TRUE
    }
  }
}
```

We can now combine the lengths we just calculated with the known species
edge lengths, this will give us the total evolutionary length and our
EDGE2 score.

``` r
# for each tip, summate ancesteral edgelengths to find EDGE2 score
for (i in 1:N_species){
  ans <- tree@edge[i,1]
  tree_data$EDGE[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
}  

tree_data$ED <- tree_data$EDGE / tree_data$pext
tree_nj <- as(tree, "phylo")
edge.res <- list(tree_data,tree)
```

To find our conservation priorities we should get a list of the highest
EDGE scores, to do this we can just arrange by descending EDGE scores.

It’s important that we also consider the other recorded data for each
species, we can see from the data that there is a spread of both
endangered species and unique species. This should allow us to draw
conclusions on current conservation priority techniques.

``` r
#print the top 10 species by highest EDGE
EDGEresults2 <- tree_data %>%
  arrange(desc(EDGE))
head(EDGEresults2, n=10)
```

    ##                       species       TBL     pext        ED       EDGE
    ## 110-10       Numenius_arquata 1.2072123 0.242500 1.2754393 0.30929403
    ## 132-133     Sterna_paradisaea 0.7721445 0.121250 0.8358817 0.10135066
    ## 107-108         Ardea_cinerea 1.2513157 0.060625 1.2842274 0.07785628
    ## 143-144      Rallus_aquaticus 1.0799354 0.060625 1.0799354 0.06547109
    ## 123-124             Alle_alle 0.3861244 0.121250 0.3994967 0.04843897
    ## 125-128     Vanellus_vanellus 0.1937066 0.242500 0.1937066 0.04697386
    ## 140-141 Streptopelia_decaocto 0.7301861 0.060625 0.7526972 0.04563227
    ## 107-5             Gavia_immer 0.6694892 0.060625 0.7024009 0.04258305
    ## 147-148       Hirundo_rustica 0.6786628 0.060625 0.6960317 0.04219692
    ## 153-16         Morus_bassanus 0.6796955 0.060625 0.6797546 0.04121012

``` r
# Find the spread of EDGE values
hist(tree_data$EDGE, breaks = 20)
```

```{image} practical_3_files/figure-gfm/unnamed-chunk-20-1.png
:align: center
:width: 600px
```

Iceland has a very skewed distribution of EDGE scores, this could be due
to a few reasons.

It’s likely due to a lack of functional variation seen in other
countries as well as a healthy overall bird population. Countries with
many more endangered species or extremely unique species will show a
better spread.

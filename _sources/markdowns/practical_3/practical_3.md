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

``` r
conservation_data$BirdLife.name = conservation_data$ScientificName
left_join(conservation_data,crosswalk, by = 'BirdLife.name')
```

    ##                      CommonName            ScientificName                             FamilyName IUCNRedListCategory BreedingEndemic
    ## 1                Rock Ptarmigan              Lagopus muta Pheasants, Partridges, Turkeys, Grouse       Least Concern           False
    ## 2                  Whooper Swan             Cygnus cygnus                    Ducks, Geese, Swans       Least Concern           False
    ## 3                   Tundra Swan        Cygnus columbianus                    Ducks, Geese, Swans       Least Concern           False
    ## 4                Barnacle Goose          Branta leucopsis                    Ducks, Geese, Swans       Least Concern           False
    ## 5                 Greylag Goose               Anser anser                    Ducks, Geese, Swans       Least Concern           False
    ## 6             Pink-footed Goose      Anser brachyrhynchus                    Ducks, Geese, Swans       Least Concern           False
    ## 7   Greater White-fronted Goose           Anser albifrons                    Ducks, Geese, Swans       Least Concern           False
    ## 8              Long-tailed Duck         Clangula hyemalis                    Ducks, Geese, Swans          Vulnerable           False
    ## 9                    King Eider     Somateria spectabilis                    Ducks, Geese, Swans       Least Concern           False
    ## 10                 Common Eider      Somateria mollissima                    Ducks, Geese, Swans     Near Threatened           False
    ## 11                Common Scoter           Melanitta nigra                    Ducks, Geese, Swans       Least Concern           False
    ## 12             Common Goldeneye        Bucephala clangula                    Ducks, Geese, Swans       Least Concern           False
    ## 13           Barrow's Goldeneye       Bucephala islandica                    Ducks, Geese, Swans       Least Concern           False
    ## 14                    Goosander          Mergus merganser                    Ducks, Geese, Swans       Least Concern           False
    ## 15       Red-breasted Merganser           Mergus serrator                    Ducks, Geese, Swans       Least Concern           False
    ## 16               Harlequin Duck Histrionicus histrionicus                    Ducks, Geese, Swans       Least Concern           False
    ## 17              Common Shelduck           Tadorna tadorna                    Ducks, Geese, Swans       Least Concern           False
    ## 18               Common Pochard             Aythya ferina                    Ducks, Geese, Swans          Vulnerable           False
    ## 19                  Tufted Duck           Aythya fuligula                    Ducks, Geese, Swans       Least Concern           False
    ## 20                Greater Scaup             Aythya marila                    Ducks, Geese, Swans       Least Concern           False
    ## 21            Northern Shoveler          Spatula clypeata                    Ducks, Geese, Swans       Least Concern           False
    ## 22                      Gadwall           Mareca strepera                    Ducks, Geese, Swans       Least Concern           False
    ## 23              Eurasian Wigeon           Mareca penelope                    Ducks, Geese, Swans       Least Concern           False
    ## 24                      Mallard        Anas platyrhynchos                    Ducks, Geese, Swans       Least Concern           False
    ## 25             Northern Pintail                Anas acuta                    Ducks, Geese, Swans       Least Concern           False
    ## 26                  Common Teal               Anas crecca                    Ducks, Geese, Swans       Least Concern           False
    ## 27                 Horned Grebe          Podiceps auritus                                 Grebes          Vulnerable           False
    ## 28                    Rock Dove             Columba livia                         Pigeons, Doves       Least Concern           False
    ## 29            Common Woodpigeon          Columba palumbus                         Pigeons, Doves       Least Concern           False
    ## 30       Eurasian Collared-dove     Streptopelia decaocto                         Pigeons, Doves       Least Concern           False
    ## 31           Western Water Rail          Rallus aquaticus               Rails, Gallinules, Coots       Least Concern           False
    ## 32                Eurasian Coot               Fulica atra               Rails, Gallinules, Coots       Least Concern           False
    ## 33            Red-throated Loon            Gavia stellata                           Loons/Divers       Least Concern           False
    ## 34                  Common Loon               Gavia immer                           Loons/Divers       Least Concern           False
    ## 35        European Storm-petrel      Hydrobates pelagicus                 Northern Storm-petrels       Least Concern           False
    ## 36         Leach's Storm-petrel     Hydrobates leucorhous                 Northern Storm-petrels          Vulnerable           False
    ## 37              Northern Fulmar        Fulmarus glacialis                   Petrels, Shearwaters       Least Concern           False
    ## 38              Manx Shearwater         Puffinus puffinus                   Petrels, Shearwaters       Least Concern           False
    ## 39                   Grey Heron             Ardea cinerea                                 Herons       Least Concern           False
    ## 40              Northern Gannet            Morus bassanus                       Gannets, Boobies       Least Concern           False
    ## 41                European Shag       Gulosus aristotelis                             Cormorants       Least Concern           False
    ## 42              Great Cormorant       Phalacrocorax carbo                             Cormorants       Least Concern           False
    ## 43       Eurasian Oystercatcher     Haematopus ostralegus                         Oystercatchers     Near Threatened           False
    ## 44                  Grey Plover      Pluvialis squatarola                                Plovers          Vulnerable           False
    ## 45       Eurasian Golden Plover       Pluvialis apricaria                                Plovers       Least Concern           False
    ## 46         Common Ringed Plover      Charadrius hiaticula                                Plovers       Least Concern           False
    ## 47             Northern Lapwing         Vanellus vanellus                                Plovers     Near Threatened           False
    ## 48                     Whimbrel         Numenius phaeopus         Sandpipers, Snipes, Phalaropes       Least Concern           False
    ## 49              Eurasian Curlew          Numenius arquata         Sandpipers, Snipes, Phalaropes     Near Threatened           False
    ## 50          Black-tailed Godwit             Limosa limosa         Sandpipers, Snipes, Phalaropes     Near Threatened           False
    ## 51              Ruddy Turnstone        Arenaria interpres         Sandpipers, Snipes, Phalaropes     Near Threatened           False
    ## 52                     Red Knot          Calidris canutus         Sandpipers, Snipes, Phalaropes     Near Threatened           False
    ## 53                   Sanderling             Calidris alba         Sandpipers, Snipes, Phalaropes       Least Concern           False
    ## 54                       Dunlin           Calidris alpina         Sandpipers, Snipes, Phalaropes     Near Threatened           False
    ## 55             Purple Sandpiper         Calidris maritima         Sandpipers, Snipes, Phalaropes       Least Concern           False
    ## 56                 Common Snipe       Gallinago gallinago         Sandpipers, Snipes, Phalaropes       Least Concern           False
    ## 57         Red-necked Phalarope        Phalaropus lobatus         Sandpipers, Snipes, Phalaropes       Least Concern           False
    ## 58                Red Phalarope     Phalaropus fulicarius         Sandpipers, Snipes, Phalaropes       Least Concern           False
    ## 59              Common Redshank            Tringa totanus         Sandpipers, Snipes, Phalaropes       Least Concern           False
    ## 60               Wood Sandpiper           Tringa glareola         Sandpipers, Snipes, Phalaropes       Least Concern           False
    ## 61                Sabine's Gull               Xema sabini                 Gulls, Terns, Skimmers       Least Concern           False
    ## 62       Black-legged Kittiwake          Rissa tridactyla                 Gulls, Terns, Skimmers          Vulnerable           False
    ## 63            Black-headed Gull          Larus ridibundus                 Gulls, Terns, Skimmers       Least Concern           False
    ## 64                     Mew Gull               Larus canus                 Gulls, Terns, Skimmers       Least Concern           False
    ## 65     Lesser Black-backed Gull              Larus fuscus                 Gulls, Terns, Skimmers       Least Concern           False
    ## 66        European Herring Gull          Larus argentatus                 Gulls, Terns, Skimmers       Least Concern           False
    ## 67                 Iceland Gull          Larus glaucoides                 Gulls, Terns, Skimmers       Least Concern           False
    ## 68                 Iceland Gull          Larus glaucoides                 Gulls, Terns, Skimmers       Least Concern           False
    ## 69                Glaucous Gull         Larus hyperboreus                 Gulls, Terns, Skimmers       Least Concern           False
    ## 70      Great Black-backed Gull             Larus marinus                 Gulls, Terns, Skimmers       Least Concern           False
    ## 71                  Arctic Tern         Sterna paradisaea                 Gulls, Terns, Skimmers       Least Concern           False
    ## 72           Long-tailed Jaeger  Stercorarius longicaudus                                  Skuas       Least Concern           False
    ## 73                Arctic Jaeger  Stercorarius parasiticus                                  Skuas       Least Concern           False
    ## 74              Pomarine Jaeger    Stercorarius pomarinus                                  Skuas       Least Concern           False
    ## 75                   Great Skua           Catharacta skua                                  Skuas       Least Concern           False
    ## 76              Atlantic Puffin        Fratercula arctica                                   Auks          Vulnerable           False
    ## 77              Black Guillemot            Cepphus grylle                                   Auks       Least Concern           False
    ## 78                    Razorbill                Alca torda                                   Auks       Least Concern           False
    ## 79                    Great Auk        Pinguinus impennis                                   Auks             Extinct           False
    ## 80                   Little Auk                 Alle alle                                   Auks       Least Concern           False
    ## 81           Thick-billed Murre               Uria lomvia                                   Auks       Least Concern           False
    ## 82                 Common Murre                Uria aalge                                   Auks       Least Concern           False
    ## 83               Long-eared Owl                 Asio otus                           Typical Owls       Least Concern           False
    ## 84              Short-eared Owl             Asio flammeus                           Typical Owls       Least Concern           False
    ## 85                    Snowy Owl           Bubo scandiacus                           Typical Owls          Vulnerable           False
    ## 86       White-tailed Sea-eagle      Haliaeetus albicilla                          Hawks, Eagles       Least Concern           False
    ## 87                       Merlin         Falco columbarius                     Falcons, Caracaras       Least Concern           False
    ## 88                    Gyrfalcon          Falco rusticolus                     Falcons, Caracaras       Least Concern           False
    ## 89                 Common Raven              Corvus corax                         Crows and jays       Least Concern           False
    ## 90        Northern House Martin          Delichon urbicum                   Swallows and martins       Least Concern           False
    ## 91                 Barn Swallow           Hirundo rustica                   Swallows and martins       Least Concern           False
    ## 92               Willow Warbler    Phylloscopus trochilus                          Leaf-warblers       Least Concern           False
    ## 93            Common Chiffchaff    Phylloscopus collybita                          Leaf-warblers       Least Concern           False
    ## 94            Eurasian Blackcap        Sylvia atricapilla                     Old World Warblers       Least Concern           False
    ## 95                Northern Wren   Troglodytes troglodytes                                  Wrens       Least Concern           False
    ## 96              Common Starling          Sturnus vulgaris                              Starlings       Least Concern           False
    ## 97                      Redwing            Turdus iliacus                               Thrushes       Least Concern           False
    ## 98           Eurasian Blackbird             Turdus merula                               Thrushes       Least Concern           False
    ## 99                    Fieldfare            Turdus pilaris                               Thrushes       Least Concern           False
    ## 100           Northern Wheatear         Oenanthe oenanthe        Old World Flycatchers and Chats       Least Concern           False
    ## 101               House Sparrow         Passer domesticus                     Old World Sparrows       Least Concern           False
    ## 102                Meadow Pipit          Anthus pratensis                    Pipits and Wagtails       Least Concern           False
    ## 103                  Rock Pipit           Anthus petrosus                    Pipits and Wagtails       Least Concern           False
    ## 104                Grey Wagtail         Motacilla cinerea                    Pipits and Wagtails       Least Concern           False
    ## 105               White Wagtail            Motacilla alba                    Pipits and Wagtails       Least Concern           False
    ## 106            Common Chaffinch         Fringilla coelebs                                Finches       Least Concern           False
    ## 107                   Brambling  Fringilla montifringilla                                Finches       Least Concern           False
    ## 108                     Redpoll          Acanthis flammea                                Finches       Least Concern           False
    ## 109                     Redpoll          Acanthis flammea                                Finches       Least Concern           False
    ## 110             Eurasian Siskin             Spinus spinus                                Finches       Least Concern           False
    ## 111            Lapland Longspur      Calcarius lapponicus                              Longspurs       Least Concern           False
    ## 112                Snow Bunting     Plectrophenax nivalis                              Longspurs       Least Concern           False
    ##                 BirdLife.name                 Jetz.name       Match.type
    ## 1                Lagopus muta              Lagopus muta     1BL to 1Jetz
    ## 2               Cygnus cygnus             Cygnus cygnus     1BL to 1Jetz
    ## 3          Cygnus columbianus        Cygnus columbianus     1BL to 1Jetz
    ## 4            Branta leucopsis          Branta leucopsis     1BL to 1Jetz
    ## 5                 Anser anser               Anser anser     1BL to 1Jetz
    ## 6        Anser brachyrhynchus      Anser brachyrhynchus     1BL to 1Jetz
    ## 7             Anser albifrons           Anser albifrons     1BL to 1Jetz
    ## 8           Clangula hyemalis         Clangula hyemalis     1BL to 1Jetz
    ## 9       Somateria spectabilis     Somateria spectabilis     1BL to 1Jetz
    ## 10       Somateria mollissima      Somateria mollissima     1BL to 1Jetz
    ## 11            Melanitta nigra           Melanitta nigra Many BL to 1Jetz
    ## 12         Bucephala clangula        Bucephala clangula     1BL to 1Jetz
    ## 13        Bucephala islandica       Bucephala islandica     1BL to 1Jetz
    ## 14           Mergus merganser          Mergus merganser     1BL to 1Jetz
    ## 15            Mergus serrator           Mergus serrator     1BL to 1Jetz
    ## 16  Histrionicus histrionicus Histrionicus histrionicus     1BL to 1Jetz
    ## 17            Tadorna tadorna           Tadorna tadorna     1BL to 1Jetz
    ## 18              Aythya ferina             Aythya ferina     1BL to 1Jetz
    ## 19            Aythya fuligula           Aythya fuligula     1BL to 1Jetz
    ## 20              Aythya marila             Aythya marila     1BL to 1Jetz
    ## 21           Spatula clypeata             Anas clypeata     1BL to 1Jetz
    ## 22            Mareca strepera             Anas strepera     1BL to 1Jetz
    ## 23            Mareca penelope             Anas penelope     1BL to 1Jetz
    ## 24         Anas platyrhynchos        Anas platyrhynchos     1BL to 1Jetz
    ## 25                 Anas acuta                Anas acuta     1BL to 1Jetz
    ## 26                Anas crecca               Anas crecca Many BL to 1Jetz
    ## 27           Podiceps auritus          Podiceps auritus     1BL to 1Jetz
    ## 28              Columba livia             Columba livia     1BL to 1Jetz
    ## 29           Columba palumbus          Columba palumbus     1BL to 1Jetz
    ## 30      Streptopelia decaocto     Streptopelia decaocto Many BL to 1Jetz
    ## 31           Rallus aquaticus          Rallus aquaticus Many BL to 1Jetz
    ## 32                Fulica atra               Fulica atra     1BL to 1Jetz
    ## 33             Gavia stellata            Gavia stellata     1BL to 1Jetz
    ## 34                Gavia immer               Gavia immer     1BL to 1Jetz
    ## 35       Hydrobates pelagicus      Hydrobates pelagicus     1BL to 1Jetz
    ## 36      Hydrobates leucorhous     Oceanodroma leucorhoa Many BL to 1Jetz
    ## 37         Fulmarus glacialis        Fulmarus glacialis     1BL to 1Jetz
    ## 38          Puffinus puffinus         Puffinus puffinus     1BL to 1Jetz
    ## 39              Ardea cinerea             Ardea cinerea     1BL to 1Jetz
    ## 40             Morus bassanus            Morus bassanus     1BL to 1Jetz
    ## 41        Gulosus aristotelis Phalacrocorax aristotelis     1BL to 1Jetz
    ## 42        Phalacrocorax carbo       Phalacrocorax carbo     1BL to 1Jetz
    ## 43      Haematopus ostralegus     Haematopus ostralegus     1BL to 1Jetz
    ## 44       Pluvialis squatarola      Pluvialis squatarola     1BL to 1Jetz
    ## 45        Pluvialis apricaria       Pluvialis apricaria     1BL to 1Jetz
    ## 46       Charadrius hiaticula      Charadrius hiaticula     1BL to 1Jetz
    ## 47          Vanellus vanellus         Vanellus vanellus     1BL to 1Jetz
    ## 48          Numenius phaeopus         Numenius phaeopus     1BL to 1Jetz
    ## 49           Numenius arquata          Numenius arquata     1BL to 1Jetz
    ## 50              Limosa limosa             Limosa limosa     1BL to 1Jetz
    ## 51         Arenaria interpres        Arenaria interpres     1BL to 1Jetz
    ## 52           Calidris canutus          Calidris canutus     1BL to 1Jetz
    ## 53              Calidris alba             Calidris alba     1BL to 1Jetz
    ## 54            Calidris alpina           Calidris alpina     1BL to 1Jetz
    ## 55          Calidris maritima         Calidris maritima     1BL to 1Jetz
    ## 56        Gallinago gallinago       Gallinago gallinago Many BL to 1Jetz
    ## 57         Phalaropus lobatus        Phalaropus lobatus     1BL to 1Jetz
    ## 58      Phalaropus fulicarius     Phalaropus fulicarius     1BL to 1Jetz
    ## 59             Tringa totanus            Tringa totanus     1BL to 1Jetz
    ## 60            Tringa glareola           Tringa glareola     1BL to 1Jetz
    ## 61                Xema sabini               Xema sabini     1BL to 1Jetz
    ## 62           Rissa tridactyla          Rissa tridactyla     1BL to 1Jetz
    ## 63           Larus ridibundus          Larus ridibundus     1BL to 1Jetz
    ## 64                Larus canus               Larus canus     1BL to 1Jetz
    ## 65               Larus fuscus              Larus fuscus     1BL to 1Jetz
    ## 66           Larus argentatus          Larus argentatus Many BL to 1Jetz
    ## 67           Larus glaucoides          Larus glaucoides 1BL to many Jetz
    ## 68           Larus glaucoides             Larus thayeri 1BL to many Jetz
    ## 69          Larus hyperboreus         Larus hyperboreus     1BL to 1Jetz
    ## 70              Larus marinus             Larus marinus     1BL to 1Jetz
    ## 71          Sterna paradisaea         Sterna paradisaea     1BL to 1Jetz
    ## 72   Stercorarius longicaudus  Stercorarius longicaudus     1BL to 1Jetz
    ## 73   Stercorarius parasiticus  Stercorarius parasiticus     1BL to 1Jetz
    ## 74     Stercorarius pomarinus    Stercorarius pomarinus     1BL to 1Jetz
    ## 75            Catharacta skua           Catharacta skua     1BL to 1Jetz
    ## 76         Fratercula arctica        Fratercula arctica     1BL to 1Jetz
    ## 77             Cepphus grylle            Cepphus grylle     1BL to 1Jetz
    ## 78                 Alca torda                Alca torda     1BL to 1Jetz
    ## 79         Pinguinus impennis                      <NA>          Extinct
    ## 80                  Alle alle                 Alle alle     1BL to 1Jetz
    ## 81                Uria lomvia               Uria lomvia     1BL to 1Jetz
    ## 82                 Uria aalge                Uria aalge     1BL to 1Jetz
    ## 83                  Asio otus                 Asio otus     1BL to 1Jetz
    ## 84              Asio flammeus             Asio flammeus     1BL to 1Jetz
    ## 85            Bubo scandiacus            Bubo scandiaca     1BL to 1Jetz
    ## 86       Haliaeetus albicilla      Haliaeetus albicilla     1BL to 1Jetz
    ## 87          Falco columbarius         Falco columbarius     1BL to 1Jetz
    ## 88           Falco rusticolus          Falco rusticolus     1BL to 1Jetz
    ## 89               Corvus corax              Corvus corax     1BL to 1Jetz
    ## 90           Delichon urbicum          Delichon urbicum Many BL to 1Jetz
    ## 91            Hirundo rustica           Hirundo rustica     1BL to 1Jetz
    ## 92     Phylloscopus trochilus    Phylloscopus trochilus     1BL to 1Jetz
    ## 93     Phylloscopus collybita    Phylloscopus collybita Many BL to 1Jetz
    ## 94         Sylvia atricapilla        Sylvia atricapilla     1BL to 1Jetz
    ## 95    Troglodytes troglodytes   Troglodytes troglodytes Many BL to 1Jetz
    ## 96           Sturnus vulgaris          Sturnus vulgaris     1BL to 1Jetz
    ## 97             Turdus iliacus            Turdus iliacus     1BL to 1Jetz
    ## 98              Turdus merula             Turdus merula Many BL to 1Jetz
    ## 99             Turdus pilaris            Turdus pilaris     1BL to 1Jetz
    ## 100         Oenanthe oenanthe         Oenanthe oenanthe Many BL to 1Jetz
    ## 101         Passer domesticus         Passer domesticus Many BL to 1Jetz
    ## 102          Anthus pratensis          Anthus pratensis     1BL to 1Jetz
    ## 103           Anthus petrosus           Anthus petrosus     1BL to 1Jetz
    ## 104         Motacilla cinerea         Motacilla cinerea     1BL to 1Jetz
    ## 105            Motacilla alba            Motacilla alba     1BL to 1Jetz
    ## 106         Fringilla coelebs         Fringilla coelebs     1BL to 1Jetz
    ## 107  Fringilla montifringilla  Fringilla montifringilla     1BL to 1Jetz
    ## 108          Acanthis flammea      Carduelis hornemanni 1BL to many Jetz
    ## 109          Acanthis flammea         Carduelis flammea 1BL to many Jetz
    ## 110             Spinus spinus          Carduelis spinus     1BL to 1Jetz
    ## 111      Calcarius lapponicus      Calcarius lapponicus     1BL to 1Jetz
    ## 112     Plectrophenax nivalis     Plectrophenax nivalis     1BL to 1Jetz

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

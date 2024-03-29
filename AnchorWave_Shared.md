### Manisha has translated the AnchorWave output to determine regions of the B73 genome that are shared or unique in each other NAM line.  

Here's what I'm going to do:  
Take my list of 175 candidate B73 NLRs (from the NAM paper) and ask if each of these falls in a shared, polymorphic, or ambiguous region \
against each NAM genome.  This will tell me something about conservation of each NLR.  

Next I will create a coordinate list of 2kb up and 2kb down from each NLR and intersect this with TE annotation to determine if these are \
share TEs.  
I will first use BEDtools merge to merge overlapping intervals (these are in effect, NLR clusters), followed by BEDtools intersect with TEs. 


<details><summary>Open to view Data</summary>
  </br>

  | V1    | V2        | V3        | V4              | V5   | NLR | KKNewList   | Distance between NLRs | Cluster(3kb) | 2kbflank1 | 2kbflank2 |
|-------|-----------|-----------|-----------------|------|-----|-------------|-----------------------|--------------|-----------|-----------|
| chr1  | 35063891  | 35068791  | Zm00001eb010990 | Gene | NLR | KK_new_list | 17783294              | NO           | 35061891  | 35070791  |
| chr1  | 52852085  | 52855107  | Zm00001eb015450 | Gene | NLR | KK_new_list | 46029646              | NO           | 52850085  | 52857107  |
| chr1  | 98884753  | 98889758  | Zm00001eb024240 | Gene | NLR | KK_new_list | 69097442              | NO           | 98882753  | 98891758  |
| chr1  | 167987200 | 167991827 | Zm00001eb030390 | Gene | NLR | KK_new_list | 16436246              | NO           | 167985200 | 167993827 |
| chr1  | 184428073 | 184439205 | Zm00001eb033040 | Gene | NLR | KK_new_list | 16603702              | NO           | 184426073 | 184441205 |
| chr1  | 201042907 | 201046494 | Zm00001eb037410 | Gene | NLR | KK_new_list | 25125093              | NO           | 201040907 | 201048494 |
| chr1  | 226171587 | 226177701 | Zm00001eb042890 | Gene | NLR | KK_new_list | 4731839               | NO           | 226169587 | 226179701 |
| chr1  | 230909540 | 230913731 | Zm00001eb044020 | Gene | NLR | KK_new_list | 66776409              | NO           | 230907540 | 230915731 |
| chr1  | 297690140 | 297701616 | Zm00001eb061810 | Gene | NLR | KK_new_list | 4192227               | NO           | 297688140 | 297703616 |
| chr1  | 301893843 | 301896543 | Zm00001eb063200 | Gene | NLR | #N/A        | 304                   | CLUSTER      | 301891843 | 301898543 |
| chr1  | 301896847 | 301904989 | Zm00001eb063210 | Gene | NLR | #N/A        | NO                    | NO           | 301894847 | 301906989 |
| chr10 | 1566941   | 1573505   | Zm00001eb405270 | Gene | NLR | KK_new_list | 68196                 | NO           | 1564941   | 1575505   |
| chr10 | 1641701   | 1648836   | Zm00001eb405290 | Gene | NLR | KK_new_list | 123082                | NO           | 1639701   | 1650836   |
| chr10 | 1771918   | 1776181   | Zm00001eb405370 | Gene | NLR | KK_new_list | 156652                | NO           | 1769918   | 1778181   |
| chr10 | 1932833   | 1934044   | Zm00001eb405380 | Gene | NLR | #N/A        | 460043                | NO           | 1930833   | 1936044   |
| chr10 | 2394087   | 2400185   | Zm00001eb405670 | Gene | NLR | REAL        | 120346                | NO           | 2392087   | 2402185   |
| chr10 | 2520531   | 2524130   | Zm00001eb405700 | Gene | NLR | KK_new_list | 180163                | NO           | 2518531   | 2526130   |
| chr10 | 2704293   | 2715986   | Zm00001eb405770 | Gene | NLR | KK_new_list | 106322                | NO           | 2702293   | 2717986   |
| chr10 | 2822308   | 2827757   | Zm00001eb405860 | Gene | NLR | KK_new_list | 5174                  | NO           | 2820308   | 2829757   |
| chr10 | 2832931   | 2844283   | Zm00001eb405870 | Gene | NLR | KK_new_list | 67247                 | NO           | 2830931   | 2846283   |
| chr10 | 2911530   | 3233043   | Zm00001eb405880 | Gene | NLR | KK_new_list | 1637                  | CLUSTER      | 2909530   | 3235043   |
| chr10 | 3234680   | 3236152   | Zm00001eb405900 | Gene | NLR | KK_new_list | 11522                 | NO           | 3232680   | 3238152   |
| chr10 | 3247674   | 3253625   | Zm00001eb405910 | Gene | NLR | KK_new_list | 64488                 | NO           | 3245674   | 3255625   |
| chr10 | 3318113   | 3319095   | Zm00001eb405920 | Gene | NLR | #N/A        | 14844                 | NO           | 3316113   | 3321095   |
| chr10 | 3333939   | 3335714   | Zm00001eb405930 | Gene | NLR | KK_new_list | 63865                 | NO           | 3331939   | 3337714   |
| chr10 | 3399579   | 3404352   | Zm00001eb405940 | Gene | NLR | KK_new_list | 37416                 | NO           | 3397579   | 3406352   |
| chr10 | 3441768   | 3447157   | Zm00001eb405960 | Gene | NLR | KK_new_list | 79174                 | NO           | 3439768   | 3449157   |
| chr10 | 3526331   | 3532500   | Zm00001eb405980 | Gene | NLR | KK_new_list | 1056140               | NO           | 3524331   | 3534500   |
| chr10 | 4588640   | 4593178   | Zm00001eb406540 | Gene | NLR | KK_new_list | 7346                  | NO           | 4586640   | 4595178   |
| chr10 | 4600524   | 4602520   | Zm00001eb406550 | Gene | NLR | #N/A        | 5138170               | NO           | 4598524   | 4604520   |
| chr10 | 9740690   | 9743607   | Zm00001eb407930 | Gene | NLR | KK_new_list | 492                   | CLUSTER      | 9738690   | 9745607   |
| chr10 | 9744099   | 9760530   | Zm00001eb407940 | Gene | NLR | KK_new_list | 18368675              | NO           | 9742099   | 9762530   |
| chr10 | 28129205  | 28134407  | Zm00001eb410900 | Gene | NLR | KK_new_list | 33080943              | NO           | 28127205  | 28136407  |
| chr10 | 61215350  | 61217574  | Zm00001eb413380 | Gene | NLR | #N/A        | 9403365               | NO           | 61213350  | 61219574  |
| chr10 | 70620939  | 70627663  | Zm00001eb414490 | Gene | NLR | KK_new_list | 14696668              | NO           | 70618939  | 70629663  |
| chr10 | 85324331  | 85329344  | Zm00001eb416830 | Gene | NLR | KK_new_list | 11610955              | NO           | 85322331  | 85331344  |
| chr10 | 96940299  | 96944201  | Zm00001eb418840 | Gene | NLR | KK_new_list | 2429330               | NO           | 96938299  | 96946201  |
| chr10 | 99373531  | 99389974  | Zm00001eb419270 | Gene | NLR | KK_new_list | 56800                 | NO           | 99371531  | 99391974  |
| chr10 | 99446774  | 99453384  | Zm00001eb419300 | Gene | NLR | #N/A        | 3870                  | NO           | 99444774  | 99455384  |
| chr10 | 99457254  | 99458927  | Zm00001eb419310 | Gene | NLR | #N/A        | 2863                  | CLUSTER      | 99455254  | 99460927  |
| chr10 | 99461790  | 99471718  | Zm00001eb419320 | Gene | NLR | KK_new_list | 328018                | NO           | 99459790  | 99473718  |
| chr10 | 99799736  | 99813696  | Zm00001eb419360 | Gene | NLR | KK_new_list | 20820295              | NO           | 99797736  | 99815696  |
| chr10 | 120633991 | 120636293 | Zm00001eb422890 | Gene | NLR | KK_new_list | 70                    | CLUSTER      | 120631991 | 120638293 |
| chr10 | 120636363 | 120636791 | Zm00001eb422900 | Gene | NLR | #N/A        | 58                    | CLUSTER      | 120634363 | 120638791 |
| chr10 | 120636849 | 120638732 | Zm00001eb422910 | Gene | NLR | #N/A        | NO                    | NO           | 120634849 | 120640732 |
| chr2  | 30418373  | 30422931  | Zm00001eb077540 | Gene | NLR | KK_new_list | 69743611              | NO           | 30416373  | 30424931  |
| chr2  | 100166542 | 100171447 | Zm00001eb087590 | Gene | NLR | KK_new_list | 16521675              | NO           | 100164542 | 100173447 |
| chr2  | 116693122 | 116696752 | Zm00001eb089490 | Gene | NLR | KK_new_list | 21254330              | NO           | 116691122 | 116698752 |
| chr2  | 137951082 | 137962100 | Zm00001eb091490 | Gene | NLR | KK_new_list | 92688                 | NO           | 137949082 | 137964100 |
| chr2  | 138054788 | 138061170 | Zm00001eb091500 | Gene | NLR | KK_new_list | 77889619              | NO           | 138052788 | 138063170 |
| chr2  | 215950789 | 215955047 | Zm00001eb108350 | Gene | NLR | KK_new_list | 4230076               | NO           | 215948789 | 215957047 |
| chr2  | 220185123 | 220197261 | Zm00001eb110140 | Gene | NLR | KK_new_list | 759705                | NO           | 220183123 | 220199261 |
| chr2  | 220956966 | 220959768 | Zm00001eb110490 | Gene | NLR | KK_new_list | 7014141               | NO           | 220954966 | 220961768 |
| chr2  | 227973909 | 227977823 | Zm00001eb112770 | Gene | NLR | KK_new_list | 3651346               | NO           | 227971909 | 227979823 |
| chr2  | 231629169 | 231633219 | Zm00001eb113900 | Gene | NLR | KK_new_list | 3682546               | NO           | 231627169 | 231635219 |
| chr2  | 235315765 | 235322044 | Zm00001eb115030 | Gene | NLR | KK_new_list | 8073                  | NO           | 235313765 | 235324044 |
| chr2  | 235330117 | 235334226 | Zm00001eb115050 | Gene | NLR | KK_new_list | 3396722               | NO           | 235328117 | 235336226 |
| chr2  | 238730948 | 238749473 | Zm00001eb116510 | Gene | NLR | KK_new_list | 2618848               | NO           | 238728948 | 238751473 |
| chr2  | 241368321 | 241373187 | Zm00001eb117700 | Gene | NLR | KK_new_list | 19244                 | NO           | 241366321 | 241375187 |
| chr2  | 241392431 | 241398053 | Zm00001eb117720 | Gene | NLR | KK_new_list | 692789                | NO           | 241390431 | 241400053 |
| chr2  | 242090842 | 242093725 | Zm00001eb118040 | Gene | NLR | KK_new_list | 2792                  | CLUSTER      | 242088842 | 242095725 |
| chr2  | 242096517 | 242119653 | Zm00001eb118050 | Gene | NLR | #N/A        | 684709                | NO           | 242094517 | 242121653 |
| chr2  | 242804362 | 242806337 | Zm00001eb118300 | Gene | NLR | #N/A        | NO                    | NO           | 242802362 | 242808337 |
| chr3  | 61902414  | 61905310  | Zm00001eb131200 | Gene | NLR | KK_new_list | 52782073              | NO           | 61900414  | 61907310  |
| chr3  | 114687383 | 114691564 | Zm00001eb134970 | Gene | NLR | KK_new_list | 567889                | NO           | 114685383 | 114693564 |
| chr3  | 115259453 | 115262468 | Zm00001eb135090 | Gene | NLR | #N/A        | 79384                 | NO           | 115257453 | 115264468 |
| chr3  | 115341852 | 115348371 | Zm00001eb135110 | Gene | NLR | KK_new_list | 182487                | NO           | 115339852 | 115350371 |
| chr3  | 115530858 | 115534811 | Zm00001eb135130 | Gene | NLR | #N/A        | 13839316              | NO           | 115528858 | 115536811 |
| chr3  | 129374127 | 129377587 | Zm00001eb136790 | Gene | NLR | KK_new_list | 4575209               | NO           | 129372127 | 129379587 |
| chr3  | 133952796 | 133954515 | Zm00001eb137530 | Gene | NLR | KK_new_list | 214575                | NO           | 133950796 | 133956515 |
| chr3  | 134169090 | 134170927 | Zm00001eb137570 | Gene | NLR | KK_new_list | 4850105               | NO           | 134167090 | 134172927 |
| chr3  | 139021032 | 139026768 | Zm00001eb138420 | Gene | NLR | KK_new_list | 54403124              | NO           | 139019032 | 139028768 |
| chr3  | 193429892 | 193432838 | Zm00001eb150750 | Gene | NLR | KK_new_list | 143097                | NO           | 193427892 | 193434838 |
| chr3  | 193575935 | 193577493 | Zm00001eb150770 | Gene | NLR | KK_new_list | -1491                 | CLUSTER      | 193573935 | 193579493 |
| chr3  | 193576002 | 193577458 | Zm00001eb150780 | Gene | NLR | #N/A        | 1258787               | NO           | 193574002 | 193579458 |
| chr3  | 194836245 | 194840589 | Zm00001eb151150 | Gene | NLR | KK_new_list | 22606342              | NO           | 194834245 | 194842589 |
| chr3  | 217446931 | 217477147 | Zm00001eb157730 | Gene | NLR | #N/A        | NO                    | NO           | 217444931 | 217479147 |
| chr4  | 1425175   | 1438104   | Zm00001eb164570 | Gene | NLR | KK_new_list | 158725                | NO           | 1423175   | 1440104   |
| chr4  | 1596829   | 1601927   | Zm00001eb164630 | Gene | NLR | KK_new_list | 566694                | NO           | 1594829   | 1603927   |
| chr4  | 2168621   | 2172051   | Zm00001eb164870 | Gene | NLR | KK_new_list | 133642                | NO           | 2166621   | 2174051   |
| chr4  | 2305693   | 2307244   | Zm00001eb164880 | Gene | NLR | KK_new_list | 4531                  | NO           | 2303693   | 2309244   |
| chr4  | 2311775   | 2314207   | Zm00001eb164890 | Gene | NLR | #N/A        | 201196                | NO           | 2309775   | 2316207   |
| chr4  | 2515403   | 2516954   | Zm00001eb164900 | Gene | NLR | KK_new_list | 82342                 | NO           | 2513403   | 2518954   |
| chr4  | 2599296   | 2600296   | Zm00001eb164910 | Gene | NLR | #N/A        | 52855                 | NO           | 2597296   | 2602296   |
| chr4  | 2653151   | 2654841   | Zm00001eb164920 | Gene | NLR | #N/A        | 68294                 | NO           | 2651151   | 2656841   |
| chr4  | 2723135   | 2724832   | Zm00001eb164930 | Gene | NLR | #N/A        | 999                   | CLUSTER      | 2721135   | 2726832   |
| chr4  | 2725831   | 2739388   | Zm00001eb164940 | Gene | NLR | KK_new_list | 495405                | NO           | 2723831   | 2741388   |
| chr4  | 3234793   | 3255162   | Zm00001eb165170 | Gene | NLR | KK_new_list | 4721                  | NO           | 3232793   | 3257162   |
| chr4  | 3259883   | 3265271   | Zm00001eb165200 | Gene | NLR | #N/A        | 13664841              | NO           | 3257883   | 3267271   |
| chr4  | 16930112  | 16940894  | Zm00001eb169030 | Gene | NLR | KK_new_list | 27384355              | NO           | 16928112  | 16942894  |
| chr4  | 44325249  | 44345549  | Zm00001eb174770 | Gene | NLR | KK_new_list | 146565946             | NO           | 44323249  | 44347549  |
| chr4  | 190911495 | 190915285 | Zm00001eb195760 | Gene | NLR | KK_new_list | 1988149               | NO           | 190909495 | 190917285 |
| chr4  | 192903434 | 192908835 | Zm00001eb196580 | Gene | NLR | KK_new_list | 2482016               | NO           | 192901434 | 192910835 |
| chr4  | 195390851 | 195393586 | Zm00001eb197290 | Gene | NLR | KK_new_list | 8839090               | NO           | 195388851 | 195395586 |
| chr4  | 204232676 | 204237107 | Zm00001eb199520 | Gene | NLR | #N/A        | 2675069               | NO           | 204230676 | 204239107 |
| chr4  | 206912176 | 206922066 | Zm00001eb200120 | Gene | NLR | KK_new_list | 845029                | NO           | 206910176 | 206924066 |
| chr4  | 207767095 | 207769977 | Zm00001eb200420 | Gene | NLR | KK_new_list | 1671688               | NO           | 207765095 | 207771977 |
| chr4  | 209441665 | 209452490 | Zm00001eb200700 | Gene | NLR | KK_new_list | 45050                 | NO           | 209439665 | 209454490 |
| chr4  | 209497540 | 209503851 | Zm00001eb200710 | Gene | NLR | KK_new_list | 181358                | NO           | 209495540 | 209505851 |
| chr4  | 209685209 | 209708969 | Zm00001eb200740 | Gene | NLR | KK_new_list | 51628                 | NO           | 209683209 | 209710969 |
| chr4  | 209760597 | 209782508 | Zm00001eb200750 | Gene | NLR | KK_new_list | 54685                 | NO           | 209758597 | 209784508 |
| chr4  | 209837193 | 209843525 | Zm00001eb200760 | Gene | NLR | KK_new_list | 10905736              | NO           | 209835193 | 209845525 |
| chr4  | 220749261 | 220759638 | Zm00001eb202350 | Gene | NLR | KK_new_list | 3877817               | NO           | 220747261 | 220761638 |
| chr4  | 224637455 | 224643701 | Zm00001eb202940 | Gene | NLR | KK_new_list | 15781441              | NO           | 224635455 | 224645701 |
| chr4  | 240425142 | 240431109 | Zm00001eb205560 | Gene | NLR | KK_new_list | NO                    | NO           | 240423142 | 240433109 |
| chr5  | 21888981  | 21895227  | Zm00001eb219900 | Gene | NLR | #N/A        | -5933                 | CLUSTER      | 21886981  | 21897227  |
| chr5  | 21889294  | 21890637  | Zm00001eb219910 | Gene | NLR | KK_new_list | 50                    | CLUSTER      | 21887294  | 21892637  |
| chr5  | 21890687  | 21891459  | Zm00001eb219920 | Gene | NLR | #N/A        | 618                   | CLUSTER      | 21888687  | 21893459  |
| chr5  | 21892077  | 21892562  | Zm00001eb219930 | Gene | NLR | #N/A        | 19471661              | NO           | 21890077  | 21894562  |
| chr5  | 41364223  | 41365909  | Zm00001eb224260 | Gene | NLR | #N/A        | 4777                  | NO           | 41362223  | 41367909  |
| chr5  | 41370686  | 41372354  | Zm00001eb224270 | Gene | NLR | KK_new_list | 15861435              | NO           | 41368686  | 41374354  |
| chr5  | 57233789  | 57238852  | Zm00001eb226690 | Gene | NLR | KK_new_list | 116569                | NO           | 57231789  | 57240852  |
| chr5  | 57355421  | 57359985  | Zm00001eb226700 | Gene | NLR | KK_new_list | 5630                  | NO           | 57353421  | 57361985  |
| chr5  | 57365615  | 57387844  | Zm00001eb226710 | Gene | NLR | KK_new_list | 216990                | NO           | 57363615  | 57389844  |
| chr5  | 57604834  | 57609095  | Zm00001eb226720 | Gene | NLR | KK_new_list | 179769                | NO           | 57602834  | 57611095  |
| chr5  | 57788864  | 57790637  | Zm00001eb226760 | Gene | NLR | KK_new_list | 536179                | NO           | 57786864  | 57792637  |
| chr5  | 58326816  | 58330060  | Zm00001eb226880 | Gene | NLR | #N/A        | 731414                | NO           | 58324816  | 58332060  |
| chr5  | 59061474  | 59069065  | Zm00001eb227070 | Gene | NLR | KK_new_list | 7647611               | NO           | 59059474  | 59071065  |
| chr5  | 66716676  | 66721872  | Zm00001eb228790 | Gene | NLR | KK_new_list | 113625292             | NO           | 66714676  | 66723872  |
| chr5  | 180347164 | 180347727 | Zm00001eb245050 | Gene | NLR | KK_new_list | 33061524              | NO           | 180345164 | 180349727 |
| chr5  | 213409251 | 213413332 | Zm00001eb253770 | Gene | NLR | KK_new_list | NO                    | NO           | 213407251 | 213415332 |
| chr6  | 11507155  | 11511452  | Zm00001eb261200 | Gene | NLR | KK_new_list | 1298451               | NO           | 11505155  | 11513452  |
| chr6  | 12809903  | 12814283  | Zm00001eb261570 | Gene | NLR | #N/A        | 266742                | NO           | 12807903  | 12816283  |
| chr6  | 13081025  | 13084717  | Zm00001eb261610 | Gene | NLR | #N/A        | 507041                | NO           | 13079025  | 13086717  |
| chr6  | 13591758  | 13592480  | Zm00001eb261630 | Gene | NLR | #N/A        | 15                    | CLUSTER      | 13589758  | 13594480  |
| chr6  | 13592495  | 13593144  | Zm00001eb261640 | Gene | NLR | #N/A        | 122075                | NO           | 13590495  | 13595144  |
| chr6  | 13715219  | 13719719  | Zm00001eb261660 | Gene | NLR | KK_new_list | 54430318              | NO           | 13713219  | 13721719  |
| chr6  | 68150037  | 68154834  | Zm00001eb268960 | Gene | NLR | #N/A        | 19998441              | NO           | 68148037  | 68156834  |
| chr6  | 88153275  | 88157186  | Zm00001eb271410 | Gene | NLR | KK_new_list | 51457594              | NO           | 88151275  | 88159186  |
| chr6  | 139614780 | 139623752 | Zm00001eb283180 | Gene | NLR | #N/A        | 20981                 | NO           | 139612780 | 139625752 |
| chr6  | 139644733 | 139655572 | Zm00001eb283200 | Gene | NLR | KK_new_list | 27405575              | NO           | 139642733 | 139657572 |
| chr6  | 167061147 | 167065592 | Zm00001eb291370 | Gene | NLR | KK_new_list | NO                    | NO           | 167059147 | 167067592 |
| chr7  | 2369632   | 2371252   | Zm00001eb298790 | Gene | NLR | #N/A        | 39509                 | NO           | 2367632   | 2373252   |
| chr7  | 2410761   | 2412587   | Zm00001eb298800 | Gene | NLR | KK_new_list | 156471                | NO           | 2408761   | 2414587   |
| chr7  | 2569058   | 2571228   | Zm00001eb298830 | Gene | NLR | KK_new_list | 8286                  | NO           | 2567058   | 2573228   |
| chr7  | 2579514   | 2580628   | Zm00001eb298840 | Gene | NLR | #N/A        | 884                   | CLUSTER      | 2577514   | 2582628   |
| chr7  | 2581512   | 2632086   | Zm00001eb298860 | Gene | NLR | #N/A        | -37156                | CLUSTER      | 2579512   | 2634086   |
| chr7  | 2594930   | 2595683   | Zm00001eb298880 | Gene | NLR | #N/A        | 3374                  | NO           | 2592930   | 2597683   |
| chr7  | 2599057   | 2600661   | Zm00001eb298890 | Gene | NLR | KK_new_list | 23490                 | NO           | 2597057   | 2602661   |
| chr7  | 2624151   | 2624891   | Zm00001eb298920 | Gene | NLR | #N/A        | 11847                 | NO           | 2622151   | 2626891   |
| chr7  | 2636738   | 2637131   | Zm00001eb298930 | Gene | NLR | #N/A        | 79659                 | NO           | 2634738   | 2639131   |
| chr7  | 2716790   | 2719122   | Zm00001eb299040 | Gene | NLR | #N/A        | 664                   | CLUSTER      | 2714790   | 2721122   |
| chr7  | 2719786   | 2720319   | Zm00001eb299050 | Gene | NLR | #N/A        | 90995                 | NO           | 2717786   | 2722319   |
| chr7  | 2811314   | 2813424   | Zm00001eb299080 | Gene | NLR | #N/A        | 6957                  | NO           | 2809314   | 2815424   |
| chr7  | 2820381   | 2823397   | Zm00001eb299090 | Gene | NLR | #N/A        | 26252                 | NO           | 2818381   | 2825397   |
| chr7  | 2849649   | 2851885   | Zm00001eb299100 | Gene | NLR | #N/A        | 108508                | NO           | 2847649   | 2853885   |
| chr7  | 2960393   | 2961094   | Zm00001eb299160 | Gene | NLR | #N/A        | 2141611               | NO           | 2958393   | 2963094   |
| chr7  | 5102705   | 5110265   | Zm00001eb299830 | Gene | NLR | KK_new_list | -5765                 | CLUSTER      | 5100705   | 5112265   |
| chr7  | 5104500   | 5105934   | Zm00001eb299840 | Gene | NLR | KK_new_list | 24097449              | NO           | 5102500   | 5107934   |
| chr7  | 29203383  | 29247396  | Zm00001eb304830 | Gene | NLR | KK_new_list | -43750                | CLUSTER      | 29201383  | 29249396  |
| chr7  | 29203646  | 29204114  | Zm00001eb304840 | Gene | NLR | #N/A        | 43478                 | NO           | 29201646  | 29206114  |
| chr7  | 29247592  | 29260760  | Zm00001eb304860 | Gene | NLR | KK_new_list | -3574                 | CLUSTER      | 29245592  | 29262760  |
| chr7  | 29257186  | 29261106  | Zm00001eb304870 | Gene | NLR | #N/A        | 288303                | NO           | 29255186  | 29263106  |
| chr7  | 29549409  | 29552744  | Zm00001eb304920 | Gene | NLR | KK_new_list | 64680164              | NO           | 29547409  | 29554744  |
| chr7  | 94232908  | 94252757  | Zm00001eb310010 | Gene | NLR | KK_new_list | 867544                | NO           | 94230908  | 94254757  |
| chr7  | 95120301  | 95157771  | Zm00001eb310060 | Gene | NLR | KK_new_list | 51095058              | NO           | 95118301  | 95159771  |
| chr7  | 146252829 | 146267516 | Zm00001eb318600 | Gene | NLR | KK_new_list | 8841320               | NO           | 146250829 | 146269516 |
| chr7  | 155108836 | 155113412 | Zm00001eb321430 | Gene | NLR | KK_new_list | 30719                 | NO           | 155106836 | 155115412 |
| chr7  | 155144131 | 155158389 | Zm00001eb321440 | Gene | NLR | KK_new_list | 4715815               | NO           | 155142131 | 155160389 |
| chr7  | 159874204 | 159899722 | Zm00001eb322130 | Gene | NLR | KK_new_list | NO                    | NO           | 159872204 | 159901722 |
| chr8  | 29814265  | 29817242  | Zm00001eb339320 | Gene | NLR | KK_new_list | 42253927              | NO           | 29812265  | 29819242  |
| chr8  | 72071169  | 72072142  | Zm00001eb343880 | Gene | NLR | #N/A        | 2                     | CLUSTER      | 72069169  | 72074142  |
| chr8  | 72072144  | 72074447  | Zm00001eb343890 | Gene | NLR | KK_new_list | 34485784              | NO           | 72070144  | 72076447  |
| chr8  | 106560231 | 106562858 | Zm00001eb349330 | Gene | NLR | KK_new_list | 47724                 | NO           | 106558231 | 106564858 |
| chr8  | 106610582 | 106614710 | Zm00001eb349360 | Gene | NLR | KK_new_list | 28518429              | NO           | 106608582 | 106616710 |
| chr8  | 135133139 | 135142568 | Zm00001eb355090 | Gene | NLR | #N/A        | -8798                 | CLUSTER      | 135131139 | 135144568 |
| chr8  | 135133770 | 135142391 | Zm00001eb355100 | Gene | NLR | KK_new_list | 1765696               | NO           | 135131770 | 135144391 |
| chr8  | 136908087 | 136909619 | Zm00001eb355630 | Gene | NLR | KK_new_list | 22666527              | NO           | 136906087 | 136911619 |
| chr8  | 159576146 | 159578040 | Zm00001eb361650 | Gene | NLR | KK_new_list | 37661                 | NO           | 159574146 | 159580040 |
| chr8  | 159615701 | 159617260 | Zm00001eb361660 | Gene | NLR | KK_new_list | 7132087               | NO           | 159613701 | 159619260 |
| chr8  | 166749347 | 166753682 | Zm00001eb363970 | Gene | NLR | KK_new_list | NO                    | NO           | 166747347 | 166755682 |
| chr9  | 2880701   | 2881264   | Zm00001eb371700 | Gene | NLR | KK_new_list | 18158608              | NO           | 2878701   | 2883264   |
| chr9  | 21039872  | 21045940  | Zm00001eb376840 | Gene | NLR | KK_new_list | 6001849               | NO           | 21037872  | 21047940  |
| chr9  | 27047789  | 27052392  | Zm00001eb378630 | Gene | NLR | KK_new_list | 93548999              | NO           | 27045789  | 27054392  |
| chr9  | 120601391 | 120602364 | Zm00001eb391090 | Gene | NLR | #N/A        | 0                     | CLUSTER      | 120599391 | 120604364 |
| chr9  | 120602364 | 120604669 | Zm00001eb391100 | Gene | NLR | KK_new_list | NO                    | NO           | 120600364 | 120606669 |


</details>

### Manisha has classified each gene into one of three categories based on the AnchorWave alignment - 1. Shared (95% of gene falls in shared block)2. Polymorphic (95% of gene falls in polymorphic block) or 3. Ambiguous (doesn't fit the other two categories)

Each B73 vs NAM comparison looks like this: \
col5 is the percent of the gene in shared block(s) \
col6 is percent of gene in a B73 insertion seq \
col7 will always be zero in the B73 vs files \
col8 is percent of gene in unalignable block \
col9 is percent missing data \
col10 is a list of the Anchorwave blocks that correspond to the region\
col11 is the classification of the gene

````
head B73_B97_gene_classification_by_full.tsv 
id_name chr start end alignable_region structural_insertion_inB73 structural_insertion_inB97 unalignable Missing_Data AW_Blocks classification
Zm00001eb000010_T001 chr1 34617 40204 1 0 0 0 0 chr1_AW_BlockID_8 shared
Zm00001eb000020_T001 chr1 41214 46762 1 0 0 0 0 chr1_AW_BlockID_8 shared
Zm00001eb000050_T001 chr1 108554 114382 0.605181880576527 0.394303363074811 0 0 0 chr1_AW_BlockID_8,chr1_AW_BlockID_9,chr1_AW_BlockID_10,chr1_AW_BlockID_11 ambiguous
Zm00001eb000060_T001 chr1 188559 189581 1 0 0 0 0 chr1_AW_BlockID_16 shared
Zm00001eb000070_T001 chr1 190192 198832 1 0 0 0 0 chr1_AW_BlockID_16 shared
Zm00001eb000080_T001 chr1 200262 203393 1 0 0 0 0 chr1_AW_BlockID_16 shared
Zm00001eb000100_T001 chr1 206619 209723 0.999355670103093 0 0 0 0 chr1_AW_BlockID_16,chr1_AW_BlockID_17,chr1_AW_BlockID_18 shared
Zm00001eb000110_T001 chr1 246422 247242 1 0 0 0 0 chr1_AW_BlockID_22 shared
Zm00001eb000120_T001 chr1 315219 315846 1 0 0 0 0 chr1_AW_BlockID_22 shared
````

I'm 'grepping' out all the NLR rows from each comparison (I'm 100% sure there's a better way to do this)
````
(base) grep -f B73_GeneIDs.txt B73_CML247_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_CML247_classification.txt
(base) grep -f B73_GeneIDs.txt B73_CML277_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_CML277_classification.txt
(base) grep -f B73_GeneIDs.txt B73_CML322_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_CML322_classification.txt
(base) grep -f B73_GeneIDs.txt B73_CML333_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_CML333_classification.txt
(base) grep -f B73_GeneIDs.txt B73_CML52_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_CML52_classification.txt
(base) grep -f B73_GeneIDs.txt B73_CML69_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_CML69_classification.txt
(base) grep -f B73_GeneIDs.txt B73_HP301_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_HP301_classification.txt
(base) grep -f B73_GeneIDs.txt B73_IL14H_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_IL14H_classification.txt
(base) grep -f B73_GeneIDs.txt B73_Il14H_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Il14H_classification.txt
(base) grep -f B73_GeneIDs.txt B73_Ki11_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Ki11_classification.txt
(base) grep -f B73_GeneIDs.txt B73_Ki3_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Ki3_classification.txt
(base) grep -f B73_GeneIDs.txt B73_Ky21_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Ky21_classification.txt
(base) grep -f B73_GeneIDs.txt B73_M162W_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_M162W_classification.txt
(base) grep -f B73_GeneIDs.txt B73_M37W_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_M37W_classification.txt
(base) grep -f B73_GeneIDs.txt B73_Mo18W_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Mo18W_classification.txt
(base) grep -f B73_GeneIDs.txt B73_Ms71_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Ms71_classification.txt
(base) grep -f B73_GeneIDs.txt B73_NC350_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_NC350_classification.txt
(base) grep -f B73_GeneIDs.txt B73_NC358_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_NC358_classification.txt
(base) grep -f B73_GeneIDs.txt B73_Oh43_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Oh43_classification.txt
(base) grep -f B73_GeneIDs.txt B73_O7B_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Oh7B_classification.txt
(base) grep -f B73_GeneIDs.txt B73_Oh7B_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Oh7B_classification.txt
(base) grep -f B73_GeneIDs.txt B73_P39_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_P39_classification.txt
(base) grep -f B73_GeneIDs.txt B73_Tx303_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Tx303_classification.txt
(base) grep -f B73_GeneIDs.txt B73_Tzi8_gene_classification_by_full.tsv > NLR_classification/B73_NLRs_Tzi8_classification.txt
````

After this I'm going into VIM and adding a column to the end with the cultivar identifier so I can cat them into a single file
`````
#here's an example of the VIM command
#:%s/$/ B97/g

#and here's how I concatenated everything
cat *.txt > B73_NLRs_All_classification.txt
`````

Here's a table showing the counts of shared, polymorphic, and ambiguous NLR genes by cultivar (B73 vs NAM) when looking at \
Exons Only -- Note that I am filtering the polymorphic NLRs from this group because the they're too different to do SV analysis

|        | ambiguous | polymorphic | shared |
| ------ | --------- | ----------- | ------ |
| B97    | 19        | 32          | 125    |
| CML103 | 26        | 32          | 118    |
| CML228 | 23        | 33          | 120    |
| CML247 | 36        | 32          | 108    |
| CML277 | 35        | 27          | 114    |
| CML322 | 38        | 29          | 109    |
| CML333 | 16        | 43          | 117    |
| CML52  | 38        | 17          | 121    |
| CML69  | 34        | 25          | 117    |
| HP301  | 23        | 42          | 111    |
| Il14H  | 35        | 24          | 117    |
| Ki11   | 28        | 32          | 116    |
| Ki3    | 34        | 27          | 115    |
| Ky21   | 30        | 27          | 119    |
| M162W  | 27        | 33          | 116    |
| M37W   | 28        | 32          | 116    |
| Mo18W  | 30        | 22          | 124    |
| Ms71   | 28        | 28          | 120    |
| NC350  | 16        | 40          | 120    |
| NC358  | 34        | 21          | 121    |
| Oh43   | 25        | 27          | 124    |
| P39    | 33        | 24          | 119    |
| Tx303  | 21        | 41          | 114    |
| Tzi8   | 27        | 38          | 111    |


And here's a comparison when you expand beyond exons and include the entire gene sequence pluse 1kb up and downstream

|        | ambiguous | polymorphic | shared |
| ------ | --------- | ----------- | ------ |
| B97    | 76        | 27          | 73     |
| CML103 | 76        | 30          | 70     |
| CML228 | 79        | 29          | 68     |
| CML247 | 96        | 27          | 53     |
| CML277 | 94        | 20          | 62     |
| CML322 | 96        | 24          | 56     |
| CML333 | 85        | 36          | 55     |
| CML52  | 97        | 13          | 66     |
| CML69  | 93        | 19          | 64     |
| HP301  | 79        | 37          | 60     |
| Il14H  | 87        | 23          | 66     |
| Ki11   | 82        | 30          | 64     |
| Ki3    | 84        | 24          | 68     |
| Ky21   | 78        | 23          | 75     |
| M162W  | 77        | 30          | 69     |
| M37W   | 78        | 30          | 68     |
| Mo18W  | 91        | 16          | 69     |
| Ms71   | 76        | 22          | 78     |
| NC350  | 81        | 33          | 62     |
| NC358  | 87        | 20          | 69     |
| Oh43   | 74        | 24          | 78     |
| P39    | 81        | 21          | 74     |
| Tx303  | 71        | 35          | 70     |
| Tzi8   | 84        | 34          | 58     |



import sys

###
###Used to extract sample name from {pegas}:haplotype function

"""Sapmle Input 1"""
#C001_IndI
#C002_VI
#C003_JapInt
#C004_TeJ
#C005_TrJ
#C006_IndInt
#C007_IndInt
#C008_IndII
#C009_IndInt
#C010_TrJ
"""Sample Input 2"""
#obtained from attributes(class haplotype)
#$index[[1]]
#[1]   1   5   9  10  11  13  15  19  20  21  23  24  25  30  31  32  36  37  39  40  41  42
#[23]  43  45  46  49  50  53  55  58  59  60  62  66  68  72  73  84  86  94  96 107 113 114
#[45] 118 121 122 125 126 127 128 130 131 135 136 141 144 146 147 150 152 154 155 156 157 159
#[67] 161 164 165 167 168 169 170 172 173 174 175 176 177 179 181 182 184 185 188 190 192 193
#[89] 201 202 203 205 206 207 208 209 210 216 217 219 222 224 228 230 231 232 235 236 237 238
#[111] 239 246 248 250 252 253 255 259 260 263 264 266 267 271 272 273 274 276 277 278 279 281
#[133] 283 285 286 287 288 289 293 294 296 297 298 301 303 304 306 307 311 312 315 317 318 320
#[155] 323 325 326 327 328 329 336 337 340 342 344 345 347 348 352 354 355 356 357 358 359 360
#[177] 363 365 366 367 372 373 374 375 378 379 380 383 384 385 386 388 389 392 393 396 399 403


def translate_array(dictionary_list, line):
    """Translate the following to sample name"""
    """[1]   3   4  12  14  16  17  18  26  28  29  34  35  48  51  52  54  56  57  63  64  65  67"""
    mylist = line.split()
    mylist.pop(0)
    result = ''
    for i in mylist:
        result +="\t" + dictionary_list[int(i)] 
    return result


sample_list = ['n']
with open(sys.argv[1]) as header:
    for line in header:
        sample_list.append(line.rstrip())

with open(sys.argv[2]) as fh:
    HaplotypeID = ''
    translated_samples = ''
    for line in fh:
        if(line.startswith('$')):
            HaplotypeID = line.rstrip()
        elif(line.startswith('[')):
            Haplotype_line = line.rstrip()
            translated_samples += translate_array(sample_list, Haplotype_line)
        elif(line.strip() == ''):
            print(HaplotypeID + translated_samples)
            HaplotypeID = ''
            translated_samples = ''


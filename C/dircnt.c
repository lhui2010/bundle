#include <stdio.h>
#include <dirent.h>

int main(int argc, char *argv[]) {
    DIR *dir;
    struct dirent *ent;
    long count = 0;

    if( argc < 2 ){
        printf("Usage:\t dircnt [directory you want to count]\n\thttps://stackoverflow.com/questions/18649547/how-to-find-the-length-of-argv-in-c\n");
        return 1;
    }

    dir = opendir(argv[1]);

    while((ent = readdir(dir)))
            ++count;

    closedir(dir);

    printf("%s contains %ld files\n", argv[1], count);

    return 0;
}

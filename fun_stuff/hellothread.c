#include <stdio.h>
#include <pthread.h>

#define NUM_THREADS 10

void *hello (void *rank) {
    printf("Hello Thread\n");
}

int main() {
    pthread_t ids[NUM_THREADS];
    for (int i=0; i < NUM_THREADS; i++) {
        pthread_create(&ids[i], NULL, hello, &i);
    }
    for (int i=0; i < NUM_THREADS; i++) {
        pthread_join(ids[i], NULL);
    }
}

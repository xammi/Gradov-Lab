#ifndef TASK_H
#define TASK_H

#include <functional>
#include <QRunnable>

typedef std::function< void() > Lambda;

class Task : public QRunnable
{
    Lambda lambda;

private:

    void run() {
        this->lambda();
    }

public:
    Task(const Lambda & lambda) : lambda(lambda) {
        this->setAutoDelete(false);
    }
};

#endif // TASK_H

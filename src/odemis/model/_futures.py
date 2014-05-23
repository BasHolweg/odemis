# -*- coding: utf-8 -*-
'''
Created on 10 Dec 2013

@author: Éric Piel

Copyright © 2013 Éric Piel, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Odemis. If not, see http://www.gnu.org/licenses/.
'''
import collections
from concurrent import futures
from concurrent.futures._base import CANCELLED, CANCELLED_AND_NOTIFIED, FINISHED, \
    PENDING, RUNNING
from concurrent.futures.thread import ThreadPoolExecutor, _WorkItem
import logging
import time


class CancellableThreadPoolExecutor(ThreadPoolExecutor):
    """
    An extended ThreadPoolExecutor that can cancel all the jobs not yet started.
    It also allows non standard Future to be created.
    """
    def __init__(self, max_workers, cls=futures.Future):
        """
        cls (class): the type of futures to create.
        """
        ThreadPoolExecutor.__init__(self, max_workers)
        self._queue = collections.deque() # thread-safe queue of futures
        # TODO: or shall we just ask the user to give its own instance of Future
        # to submit()
        self._cls = cls

    def submit(self, fn, *args, **kwargs):
        logging.debug("queuing action %s with arguments %s", fn, args)
        with self._shutdown_lock:
            if self._shutdown:
                raise RuntimeError('cannot schedule new futures after shutdown')

            f = self._cls()
            w = _WorkItem(f, fn, args, kwargs)

            self._work_queue.put(w)
            self._adjust_thread_count()

            # add to the queue and track the task
            self._queue.append(f)
            f.add_done_callback(self._on_done)
        return f
    submit.__doc__ = ThreadPoolExecutor.submit.__doc__

    def _on_done(self, future):
        # task is over
        try:
            self._queue.remove(future)
        except ValueError:
            # can happen if it was cancelled
            pass

    def cancel(self):
        """
        Cancels all the tasks still in the work queue, if they can be cancelled
        Returns when all the tasks have been cancelled or are done.
        """
        logging.debug("Cancelling all the %d futures in queue", len(self._queue))
        uncancellables = []
        # cancel one task at a time until there is nothing in the queue
        while True:
            try:
                # Start with the last one added as it's the most likely to be cancellable
                f = self._queue.pop()
            except IndexError:
                break
            logging.debug("Cancelling %s", f)
            if not f.cancel():
                uncancellables.append(f)

        # wait for the non cancellable tasks to finish
        for f in uncancellables:
            try:
                f.result()
            except:
                # the task raised an exception => we don't care
                pass

class InstantaneousFuture(futures.Future):
    """
    This is a simple class which follows the Future interface and represents a
    call already finished when returning.
    """
    def __init__(self, result=None, exception=None):
        self._result = result
        self._exception = exception

    def cancel(self):
        return False

    def cancelled(self):
        return False

    def running(self):
        return False

    def done(self):
        return True

    def result(self, timeout=None):
        if self._exception:
            raise self._exception
        return self._result

    def exception(self, timeout=None):
        return self._exception

    def add_done_callback(self, fn):
        fn(self)


class CancellableFuture(futures.Future):
    """
    set task_canceller to a callable to allow cancelling a running task
    """

    def __init__(self):
        futures.Future.__init__(self)

        # Callable that takes the future as argument and returns True if the
        # cancellation was successful (and False otherwise).
        # As long as it's None, the future cannot be cancelled while running
        self.task_canceller = None

    def cancel(self):
        """Cancel the future if possible.

        Returns True if the future was cancelled, False otherwise. A future
        cannot be cancelled if it has already completed.
        """
        # different implementation because we _can_ cancel a running task, by
        # calling a special function
        with self._condition:
            if self._state == FINISHED:
                return False

            if self._state in [CANCELLED, CANCELLED_AND_NOTIFIED]:
                return True

            if self._state == RUNNING:
                canceller = self.task_canceller
                if not (canceller and canceller(self)):
                    return False

            self._state = CANCELLED
            self._condition.notify_all()

        self._invoke_callbacks()
        return True

class ProgressiveFuture(CancellableFuture):
    """
    Allows to track the current progress of the task by getting the (expected)
    start and end time.
    """

    def __init__(self, start=None, end=None):
        """
        start (float): start time
        end (float): end time
        """
        CancellableFuture.__init__(self)
        self._upd_callbacks = []

        # just a bit ahead of time to say it's not starting now
        self._start_time = start or (time.time() + 0.1)
        self._end_time = end or (self._start_time + 0.1)
        self.add_done_callback(self.__on_done)

    def __on_done(self, future):
        """
        Called when the future is over to report the update one last time
        => the last update callbacks will be called _before_ the done callbacks
        """
        with self._condition:
            self._end_time = time.time()
        self._invoke_upd_callbacks()

    # TODO: Add support to get the current status? as a standard call?

    def _report_update(self, fn):
        # Why the 'with'? Is there some cleanup needed?
        with self._condition:
            now = time.time()
            if self._state in [CANCELLED, CANCELLED_AND_NOTIFIED, FINISHED]:
                past = self._end_time - self._start_time
                left = 0
            elif self._state == PENDING:
                past = now - self._start_time
                # ensure we state it's not yet started
                if past >= 0:
                    past = -1e-9
                # ensure past + left == duration
                left = (self._end_time - self._start_time) - past
            else: # running
                past = now - self._start_time
                left = self._end_time - now
                if left < 0:
                    logging.debug("reporting progress on task which should have "
                                  "finished already %f s ago", -left)
                    left = 0
        try:
            # TODO: better use absolute values: start/now/end ? or current ratio/start/end?
            fn(self, past, left)
        except Exception:
            logging.exception('exception calling callback for %r', self)

    def _invoke_upd_callbacks(self):
        for callback in self._upd_callbacks:
            self._report_update(callback)

    def set_start_time(self, val):
        """
        Update the start time of the task. To be used by executors only.

        val (float): time at which the task started (or will be starting)
        """
        with self._condition:
            self._start_time = val
        self._invoke_upd_callbacks()

    def set_end_time(self, val):
        """
        Update the end time of the task. To be used by executors only.

        val (float): time at which the task ended (or will be ending)
        """
        with self._condition:
            self._end_time = val
        self._invoke_upd_callbacks()

    def add_update_callback(self, fn):
        """
        Adds a callback that will receive progress updates whenever a new one is
          available. The callback receives 2 floats: past and left.
          "past" is the number of seconds elapsed since the beginning of the
          task, and "left" is the estimated number of seconds until the end of the
          task. If the task is not yet started, past can be negative, indicating
          the estimated time before the task starts. If the task is finished (or
          cancelled) the time left is 0 and the time past is the duration of the
          task. The callback is always called at least once, when the task is
          finished.
        fn (callable: (Future, float, float) -> None): the callback, that will
          be called with this future as argument and the past and left information.
        """
        with self._condition:
            if self._state not in [CANCELLED, FINISHED]:
                self._upd_callbacks.append(fn)

        # Immediately report the current known information (even if finished)
        self._report_update(fn)

    def set_running_or_notify_cancel(self):
        """
        Returns:
            False if the Future was cancelled, True otherwise.
        """
        running = futures.Future.set_running_or_notify_cancel(self)
        if running:
            now = time.time()
            with self._condition:
                self._start_time = now
            self._invoke_upd_callbacks()

        return running



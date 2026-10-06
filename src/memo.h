/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_MEMO_H_
#define CARTA_ZARR_SRC_MEMO_H_

#include "carta-zarr/error.h"
#include "carta-zarr/result.h"

#include <list>
#include <map>
#include <mutex>
#include <utility>

namespace carta::zarr::internal {

// Whether a remembered value is to be computed again when it is next asked for: a read that could not
// be made, rather than an answer about the store. An I/O error may not happen the next time, and a
// cancelled read was stopped by its caller, not by what it read; neither says anything about the
// store that the next caller should be held to. Every other value, failures included, is an answer.
template <typename Value>
bool AskAgain(const Value& /*value*/) {
    return false;
}
template <typename T>
bool AskAgain(const Result<T>& value) {
    return !value && (value.error().code == ErrorCode::io_error || value.error().code == ErrorCode::cancelled);
}

/**
 * A table of values computed on first use and remembered thereafter.
 *
 * The lock is held across the computation, so a value is computed exactly once no matter how many
 * threads ask for it at the same time. That also means a computation must not reach back into the
 * same table: a caller whose computation consults another Memo has to order those tables and keep
 * every caller to the same order.
 *
 * A failed computation is remembered like any other value. That is deliberate for a read-only view
 * of a store: a node that was missing stays missing for the life of the view, so what a caller
 * observes does not change under it. A read that could not be made is the exception (AskAgain): it
 * is no answer about the store, so it is computed again when next asked for, until there is one.
 *
 * Values are handed back by reference. An entry is written once and never erased or replaced -- one
 * computed again is added after the one it supersedes, which stays where it was -- and neither a
 * std::map nor a std::list moves the values it already holds, so a reference stays valid for as long
 * as the table does -- which is what makes a table of parsed JSON documents worth
 * having at all. A caller that stores one keeps the table alive for at least as long.
 */
template <typename Key, typename Value>
class Memo {
public:
    Memo() = default;
    Memo(const Memo&) = delete;
    Memo& operator=(const Memo&) = delete;
    Memo(Memo&&) = delete;
    Memo& operator=(Memo&&) = delete;
    ~Memo() = default;

    template <typename Compute>
    const Value& GetOrCompute(const Key& key, Compute compute) const {
        std::scoped_lock const lock(_mutex);
        auto& answers = _entries[key];
        if (answers.empty() || AskAgain(answers.back())) {
            answers.push_back(compute());
        }
        return answers.back();
    }

    // Put a value in that was not computed here. It is for a table whose entries are already in
    // hand -- the node metadata a store's consolidated copy arrives holding -- so that they are one
    // value looked up one way rather than a second copy consulted first. An entry that is already
    // there wins, because it is the one callers may be holding a reference to.
    void Insert(Key key, Value value) const {
        std::scoped_lock const lock(_mutex);
        auto& answers = _entries[std::move(key)];
        if (answers.empty()) {
            answers.push_back(std::move(value));
        }
    }

private:
    mutable std::mutex _mutex;
    // Each key's answers in the order they were computed, the last the current one.
    mutable std::map<Key, std::list<Value>> _entries;
};

/**
 * One value computed on first use and remembered thereafter: Memo with nothing to key on.
 *
 * Shares Memo's contract, including that a failed computation is remembered unless AskAgain says it
 * is no answer.
 */
template <typename Value>
class Lazy {
public:
    Lazy() = default;
    Lazy(const Lazy&) = delete;
    Lazy& operator=(const Lazy&) = delete;
    Lazy(Lazy&&) = delete;
    Lazy& operator=(Lazy&&) = delete;
    ~Lazy() = default;

    template <typename Compute>
    const Value& GetOrCompute(Compute compute) const {
        std::scoped_lock const lock(_mutex);
        if (_answers.empty() || AskAgain(_answers.back())) {
            _answers.push_back(compute());
        }
        return _answers.back();
    }

private:
    mutable std::mutex _mutex;
    mutable std::list<Value> _answers;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_MEMO_H_
